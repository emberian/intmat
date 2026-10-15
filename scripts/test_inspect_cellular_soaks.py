import json
from pathlib import Path
import tempfile
import unittest

from inspect_cellular_soaks import checkpoint_pointer, metric_chunks


class ReadOnlyInspection(unittest.TestCase):
    def test_completed_pointer_and_path_boundary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'state.jls').write_bytes(b'opaque Julia data')
            (root/'checkpoint.json').write_text(json.dumps(dict(file='state.jls', step=17)))
            result = checkpoint_pointer(root)
            self.assertTrue(result['stable'])
            self.assertEqual(result['pointer']['step'], 17)
            self.assertEqual(result['bytes'], 17)
            self.assertEqual((root/'state.jls').read_bytes(), b'opaque Julia data')
            (root/'checkpoint.json').write_text(json.dumps(dict(file='../state.jls', step=17)))
            with self.assertRaises(ValueError):
                checkpoint_pointer(root)

    def test_incomplete_append_is_distinct_from_corrupt_complete_row(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            row = dict(step=9, metrics=dict(population=12, bonds=3, huge_internal_array=[1, 2]))
            original = json.dumps(row)+'\ninvalid\n{"in_progress":'
            (root/'metrics-one.jsonl').write_text(original)
            result = metric_chunks(root)
            self.assertEqual(result['complete_rows'], 1)
            self.assertEqual(result['invalid_complete_rows'], 1)
            self.assertTrue(result['partial_tail'])
            self.assertEqual(result['recent'][0]['metrics'], dict(population=12, bonds=3))
            self.assertEqual((root/'metrics-one.jsonl').read_text(), original)


if __name__ == '__main__':
    unittest.main()
