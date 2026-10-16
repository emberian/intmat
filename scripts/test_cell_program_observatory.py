import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('builder', Path(__file__).with_name('build_cell_program_observatory.py'))
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)

class ProgramBuilderTests(unittest.TestCase):
    def fixture(self):
        return {'schema': builder.SCHEMA, 'metadata': {'scope': 'Explicitly synthetic decoding fixture'}, 'frames': [
            {'time': .02, 'step': 1, 'cells': [{'id': 1, 'sensed_step': 0, 'raw': [.2, .3], 'raw_next': [.3, .4], 'prediction_for_step': 10}]},
            {'time': .22, 'step': 11, 'cells': [{'id': 1, 'sensed_step': 10, 'raw': [.31, .4], 'raw_next': [.4, .4], 'prediction_for_step': 20}]}]}

    def test_completed_clock_and_sensor_clock_are_preserved(self):
        document = builder.validate(self.fixture())
        self.assertEqual(document['frames'][1]['step'], 11)
        self.assertEqual(document['frames'][1]['cells'][0]['sensed_step'], 10)

    def test_invalid_identity_and_shape_are_rejected(self):
        d = self.fixture();d['frames'][0]['cells'] *= 2
        with self.assertRaises(ValueError): builder.validate(d)
        d = self.fixture();d['frames'][0]['cells'][0]['raw'] = [0]
        with self.assertRaises(ValueError): builder.validate(d)

    def test_builder_preserves_absence_and_checks_byte_bound(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp);source = root/'actual.json';d = self.fixture();d['frames'][1]['cells'] = []
            source.write_text(json.dumps(d));result = builder.build([source], root/'page')
            self.assertEqual(result['captures'], 2)
            html = (root/'page/index.html').read_text()
            self.assertIn('"cells":[]', html)
            with self.assertRaises(ValueError): builder.build([source], root/'too-small', 10)

    def test_private_source_paths_are_not_published(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp);source = root/'actual.json';d = self.fixture();d['metadata']['private'] = '/Users/example/transcript'
            source.write_text(json.dumps(d))
            with self.assertRaises(ValueError): builder.build([source], root/'page')

if __name__ == '__main__': unittest.main()
