"""Exact-byte transport, complete catalog rewriting and privacy-before-compression."""
import gzip,hashlib,json,pathlib,tempfile,unittest
from static_recording_transport import RECORDINGS,stage_recordings
class TransportTests(unittest.TestCase):
    def fixture(self,root):
        originals={}
        for name in RECORDINGS:
            path=root/name;path.parent.mkdir(parents=True,exist_ok=True)
            data=[dict(time=.12345678901234567,cells=[dict(id=17,state=[-1e-30,1.2345678901234])],label='μ')]
            if name.startswith('cellular/'):data=dict(frames=data,grid=dict(extent=[2,3]))
            raw=json.dumps(data,ensure_ascii=False).encode();path.write_bytes(raw);originals[name]=raw
        for name in ['arena','metabolic']:(root/name/'arena.html').write_text("<script>async function start(){let response=await fetch('replay.json');}</script>")
        (root/'cellular/index.html').write_text("<script>async function load(){let r=await fetch($('dataset').value);}$('dataset').value='mosaic-early.json';</script>")
        (root/'cellular/catalog.json').write_text(json.dumps([dict(file=pathlib.PurePosixPath(p).name) for p in RECORDINGS if p.startswith('cellular/')]))
        return originals
    def test_exact_bytes_states_hashes_and_references(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=pathlib.Path(temporary);originals=self.fixture(root);manifest=stage_recordings(root)
            self.assertEqual(len(manifest['recordings']),8)
            for row in manifest['recordings']:
                raw=gzip.decompress((root/row['compressed_path']).read_bytes());expected=originals[row['original_path']]
                self.assertEqual(raw,expected);self.assertEqual(json.loads(raw),json.loads(expected));self.assertEqual(hashlib.sha256(raw).hexdigest(),row['sha256']);self.assertEqual(len(raw),row['original_bytes']);self.assertFalse((root/row['original_path']).exists())
            for name in ('arena/arena.html','metabolic/arena.html','cellular/index.html'):
                text=(root/name).read_text();self.assertIn('recordingFetch',text);self.assertIn('../recording-transport.js',text)
            self.assertTrue(all(row['file'].endswith('.json.gz') for row in json.loads((root/'cellular/catalog.json').read_text())))
    def test_private_reference_rejected_before_any_mutation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=pathlib.Path(temporary);self.fixture(root);path=root/RECORDINGS[0];path.write_text('[{"source":"/Users/person/private"}]')
            with self.assertRaisesRegex(ValueError,'Private/local'):stage_recordings(root)
            self.assertTrue(path.exists());self.assertFalse(list(root.rglob('*.gz')))
    def test_unrecognized_loader_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=pathlib.Path(temporary);self.fixture(root);(root/'arena/arena.html').write_text('<script>changed()</script>')
            with self.assertRaisesRegex(ValueError,'Unexpected recording loader'):stage_recordings(root)
            self.assertTrue(all((root/p).exists() for p in RECORDINGS))
if __name__=='__main__':unittest.main()
