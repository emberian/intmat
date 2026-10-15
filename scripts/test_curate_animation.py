import hashlib
import importlib.util
import json
import pathlib
import shutil
import tempfile
import unittest

from curate_animation import curate, validate_sequence


@unittest.skipUnless(importlib.util.find_spec('PIL'), 'Pillow required for media curation')
class AnimationPublication(unittest.TestCase):
    def setUp(self):
        from PIL import Image
        self.temp = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self.temp.name)
        self.manifest = self.root / 'sequence.json'
        self.sequence = dict(id='fixture', title='Fixture', description='Synthetic encoder test only',
                             fps=4, source={'native_shape': [2, 2, 1]}, frames=[])
        for i in range(3):
            path = self.root / f'{i}.png'
            Image.new('RGB', (64, 48), (i * 90, 30, 20)).save(path)
            self.sequence['frames'].append(dict(file=path.name, time=i * .25,
                                                sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        self.approve()

    def tearDown(self):
        self.temp.cleanup()

    def approve(self):
        self.manifest.write_text(json.dumps(self.sequence))
        (self.root / 'art-reviewed.json').write_text(json.dumps(dict(status='approved',
            sequence_sha256=hashlib.sha256(self.manifest.read_bytes()).hexdigest())))

    def test_changed_frame_and_manifest_rejected(self):
        (self.root / '0.png').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'changed after review'):
            validate_sequence(self.manifest)
        self.manifest.write_text(self.manifest.read_text() + ' ')
        with self.assertRaisesRegex(ValueError, 'matching completed visual review'):
            validate_sequence(self.manifest)

    def test_unordered_time_and_escaped_path_rejected(self):
        self.sequence['frames'][1]['time'] = 0
        self.approve()
        with self.assertRaisesRegex(ValueError, 'strictly increasing'):
            validate_sequence(self.manifest)
        self.sequence['frames'][1]['time'] = .25
        self.sequence['frames'][0]['file'] = '../elsewhere.png'
        self.approve()
        with self.assertRaisesRegex(ValueError, 'inside the sequence'):
            validate_sequence(self.manifest)

    @unittest.skipUnless(shutil.which('ffmpeg') and shutil.which('ffprobe'), 'ffmpeg required')
    def test_video_roundtrip_and_budget_failure_are_atomic(self):
        self.sequence['poster_frame_index'] = 0
        self.approve()
        destination = self.root / 'public'
        with self.assertRaisesRegex(ValueError, 'exceeds publication budget'):
            curate(self.manifest, destination, budget_bytes=1)
        self.assertFalse((destination / 'fixture').exists())
        target = curate(self.manifest, destination)
        data = json.loads((target / 'film.json').read_text())
        self.assertEqual(len(data['frames']), 3)
        self.assertEqual(data['playback']['duration_seconds'], .75)
        self.assertEqual([f['simulation_time'] for f in data['frames']], [0, .25, .5])
        self.assertEqual(data['playback']['width'], 64)
        self.assertEqual(data['poster_frame_index'], 0)
        from PIL import Image
        with Image.open(target / 'poster.webp') as poster:
            self.assertLess(poster.getpixel((32, 24))[0], 10)
        with self.assertRaisesRegex(ValueError, 'Publication exists'):
            curate(self.manifest, destination)


if __name__ == '__main__':
    unittest.main()
