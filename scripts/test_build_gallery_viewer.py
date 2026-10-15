"""Exercise optional viewer metadata through a real gallery build."""
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

class GalleryViewerTest(unittest.TestCase):
    def test_optional_viewers_are_local_and_labels_are_explicit(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            def specimen(name, viewer=None):
                directory=root/name;directory.mkdir()
                metadata={'title':name,'family':'ecology'}
                if viewer is not None:metadata.update(viewer=viewer,viewer_label='Open recorded motion',viewer_context='Recorded motion + cell inspection')
                (directory/'provenance.json').write_text(json.dumps(metadata))
                (directory/'preview.png').write_bytes(b'test-only image placeholder')
                return directory
            custom=specimen('custom','index.html');(custom/'index.html').write_text('<title>Recorded motion</title>')
            legacy=specimen('legacy');(legacy/'viewer.html').write_text('<title>Volume</title>')
            nested=specimen('nested','motion/index.html');(nested/'motion').mkdir();(nested/'motion/index.html').write_text('<title>Nested motion</title>')
            outside=root/'outside.html';outside.write_text('<title>Outside specimen</title>')
            invalid=['https://example.com/index.html','javascript:alert(1)','../outside.html','/outside.html','%2e%2e/outside.html','motion\\index.html']
            for i,viewer in enumerate(invalid):specimen('invalid-'+str(i),viewer)
            linked=specimen('linked','index.html');(linked/'index.html').symlink_to(outside)
            script=Path(__file__).with_name('build_gallery.py')
            subprocess.run([sys.executable,str(script),str(root)],check=True,capture_output=True,text=True)
            match=re.search(r'<script id="specimens" type="application/json">(.*?)</script>',(root/'index.html').read_text(),re.S)
            items={item['id']:item for item in json.loads(match[1])}
            self.assertEqual(items['custom']['viewer'],'custom/index.html')
            self.assertEqual(items['custom']['viewer_label'],'Open recorded motion')
            self.assertEqual(items['custom']['viewer_context'],'Recorded motion + cell inspection')
            self.assertEqual(items['legacy']['viewer'],'legacy/viewer.html')
            self.assertEqual(items['legacy']['viewer_label'],'Open interactive volume')
            self.assertEqual(items['nested']['viewer'],'nested/motion/index.html')
            for name in ['linked']+['invalid-'+str(i) for i in range(len(invalid))]:
                self.assertIsNone(items[name]['viewer'],name)

if __name__=='__main__':unittest.main()
