"""Transport regression: preserve every quantized byte and custom viewer edit."""
import base64,json,unittest
from viewer_storage import RECORDING,compress_html,unpack

class StorageTests(unittest.TestCase):
    def test_lossless_deterministic_idempotent_and_custom_source_preserved(self):
        fields={role:[base64.b64encode(bytes(range(256))*8).decode(),base64.b64encode(bytes(2048)).decode()] for role in ['deposit','activity','nutrient']}
        data=dict(run=dict(shape=[8,16,16],frames=[dict(time=0),dict(time=2)]),scales=dict(deposit=[-1,1]),fields=fields)
        source='<script id="recording" type="application/json">'+json.dumps(data)+'</script><script>const data=JSON.parse(document.getElementById(\'recording\').textContent),gl={};if(!gl){};/* custom contour shader */</script>'
        compressed=compress_html(source)
        result=json.loads(RECORDING.search(compressed)[2])
        self.assertEqual(unpack(result),fields)
        for role in fields:
            for old,new in zip(fields[role],unpack(result)[role]):self.assertEqual(base64.b64decode(old),base64.b64decode(new))
        self.assertEqual(result['run'],data['run']);self.assertEqual(result['scales'],data['scales'])
        self.assertIn('/* custom contour shader */',compressed)
        self.assertEqual(compressed,compress_html(source));self.assertEqual(compressed,compress_html(compressed))
        self.assertLess(len(compressed),len(source))
    def test_refresh_sanitizes_and_compresses_existing_recording(self):
        import pathlib,tempfile
        import build_pages
        fields=dict(deposit=[base64.b64encode(bytes(128)).decode()])
        data=dict(run=dict(private_path='/Users/example/private'),fields=fields)
        source='<script id="recording" type="application/json">'+json.dumps(data)+'</script><script>const data=JSON.parse(document.getElementById(\'recording\').textContent),gl={};if(!gl){}</script>'
        previous=build_pages.REPO,build_pages.SOURCE,build_pages.MANIFEST
        try:
            with tempfile.TemporaryDirectory() as temporary:
                root=pathlib.Path(temporary);origin=root/'outputs'/'sample';origin.mkdir(parents=True)
                (origin/'viewer.html').write_text(source)
                manifest=root/'manifest.json';manifest.write_text(json.dumps([dict(source='sample',id='sample',files=['viewer.html'])]))
                build_pages.REPO=root;build_pages.SOURCE=root/'collection';build_pages.MANIFEST=manifest
                build_pages.refresh()
                result=json.loads(RECORDING.search((root/'collection'/'sample'/'viewer.html').read_text())[2])
                self.assertEqual(unpack(result),fields)
                self.assertEqual(result['run']['private_path'],'private')
                self.assertIn('fields_gzip',result)
        finally:
            build_pages.REPO,build_pages.SOURCE,build_pages.MANIFEST=previous
    def test_unrecognized_bootstrap_rejected(self):
        with self.assertRaises(ValueError):compress_html('<script id="recording" type="application/json">{"fields":{}}</script>')
if __name__=='__main__':unittest.main()
