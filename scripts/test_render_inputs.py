import json
import os
from pathlib import Path
import tempfile
import unittest

from render_inputs import source_signature, source_revision


class RenderInputsTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.job = {'renderer': 'interfacial-choir', 'volume': 'deliberately-unused',
                    'args': ['--steps', '10']}
        self.frames = [{'step': step, 'directory': f'frames/{step}'} for step in (0, 10)]
        (self.root / 'run.json').write_text(json.dumps({'frames': self.frames}))
        for frame in self.frames:
            directory = self.root / frame['directory']
            directory.mkdir(parents=True)
            for name in ('volume.json', 'solid.f32', 'activity.f32', 'pressure.f32'):
                (directory / name).write_bytes(b'1234')

    def test_selected_field_bytes_invalidate_even_with_preserved_timestamp(self):
        before = source_signature(self.root, self.job)
        revision = source_revision(self.root, self.job)
        path = self.root / 'frames/10/solid.f32'
        stamp = path.stat()
        path.write_bytes(b'5678')
        os.utime(path, ns=(stamp.st_atime_ns, stamp.st_mtime_ns))
        self.assertNotEqual(before, source_signature(self.root, self.job))
        self.assertNotEqual(revision, source_revision(self.root, self.job))

    def test_run_frame_mapping_is_an_input(self):
        before = source_signature(self.root, self.job)
        self.frames[1]['directory'] = 'frames/0'
        (self.root / 'run.json').write_text(json.dumps({'frames': self.frames}))
        after = source_signature(self.root, self.job)
        self.assertNotEqual(before, after)
        self.assertIn('frames/0/solid.f32', [entry[0] for entry in after])
        self.assertNotIn('frames/10/solid.f32', [entry[0] for entry in after])

    def test_unselected_raw_frame_does_not_invalidate(self):
        before = source_signature(self.root, self.job)
        (self.root / 'frames/0/pressure.f32').write_bytes(b'changed')
        self.assertEqual(before, source_signature(self.root, self.job))

    def test_per_frame_shape_manifest_is_an_input(self):
        before = source_signature(self.root, self.job)
        (self.root / 'frames/10/volume.json').write_bytes(b'new shape')
        self.assertNotEqual(before, source_signature(self.root, self.job))

    def test_equals_syntax_and_last_steps_option(self):
        self.job['args'] = ['--steps', '0', '--steps=10']
        paths = [entry[0] for entry in source_signature(self.root, self.job)]
        self.assertIn('frames/10/solid.f32', paths)
        self.assertNotIn('frames/0/solid.f32', paths)


if __name__ == '__main__':
    unittest.main()
