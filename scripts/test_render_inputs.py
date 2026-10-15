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


    def test_electrothermal_compressed_fields_and_graph_are_inputs(self):
        self.job['renderer'] = 'electrothermal-archive'
        selected = self.root / 'frames/10'
        (selected / 'geometry.json').write_text('{"edges": [[1,2]]}')
        (selected / 'temperature.f32.gz').write_bytes(b'compressed thermal field')
        signature = source_signature(self.root, self.job)
        paths = [entry[0] for entry in signature]
        self.assertIn('frames/10/geometry.json', paths)
        self.assertIn('frames/10/temperature.f32.gz', paths)
        (selected / 'temperature.f32.gz').write_bytes(b'changed thermal field')
        self.assertNotEqual(signature, source_signature(self.root, self.job))
        signature = source_signature(self.root, self.job)
        (selected / 'geometry.json').write_text('{"edges": [[2,3]]}')
        self.assertNotEqual(signature, source_signature(self.root, self.job))

    def test_electrothermal_prefers_raw_if_both_encodings_exist(self):
        self.job['renderer'] = 'electrothermal-archive'
        selected = self.root / 'frames/10'
        (selected / 'geometry.json').write_text('{}')
        (selected / 'temperature.f32').write_bytes(b'raw thermal field')
        (selected / 'temperature.f32.gz').write_bytes(b'unused compressed field')
        before = source_signature(self.root, self.job)
        (selected / 'temperature.f32.gz').write_bytes(b'changed unused compression')
        self.assertEqual(before, source_signature(self.root, self.job))


    def test_repair_weave_uses_selected_graph_without_volume(self):
        self.job['renderer'] = 'repair-weave'
        path = self.root / 'frames/10/geometry.json'
        path.write_text('{"damage": [0,1]}')
        before = source_signature(self.root, self.job)
        self.assertEqual([entry[0] for entry in before], ['run.json', 'frames/10/geometry.json'])
        (self.root / 'frames/10/solid.f32').write_bytes(b'unused')
        self.assertEqual(before, source_signature(self.root, self.job))
        path.write_text('{"damage": [1,1]}')
        self.assertNotEqual(before, source_signature(self.root, self.job))


    def test_constructive_volume_tracks_geometry_and_consumed_fields(self):
        self.job['renderer'] = 'constructive-volume'
        selected = self.root / 'frames/10'
        for name in ('geometry.json', 'retained_a.f32', 'retained_b.f32', 'permeability.f32'):
            (selected / name).write_bytes(b'original')
        before = source_signature(self.root, self.job)
        (selected / 'enzyme_a.f32').write_bytes(b'unused by optical mapping')
        self.assertEqual(before, source_signature(self.root, self.job))
        for name in ('geometry.json', 'retained_a.f32', 'retained_b.f32', 'permeability.f32'):
            before = source_signature(self.root, self.job)
            (selected / name).write_bytes(b'changed')
            self.assertNotEqual(before, source_signature(self.root, self.job))


if __name__ == '__main__':
    unittest.main()
