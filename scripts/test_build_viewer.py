#!/usr/bin/env python3
"""Small exact-byte fixtures for the legacy and streamed viewer adapters."""
import array,base64,json,pathlib,re,subprocess,sys,tempfile,unittest
from viewer_storage import unpack
SCRIPT=pathlib.Path(__file__).with_name('build_viewer.py')
class AdapterTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=pathlib.Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()
    def write(self,frame,field,values):
        v=array.array('f',values)
        if sys.byteorder!='little':v.byteswap()
        (self.root/f'{frame}-{field}.f32').write_bytes(v.tobytes())
    def build(self,*args,success=True):
        proc=subprocess.run([sys.executable,str(SCRIPT),str(self.root),*args],capture_output=True,text=True)
        if not success:self.assertNotEqual(proc.returncode,0);return proc
        self.assertEqual(proc.returncode,0,proc.stderr)
        text=(self.root/'viewer.html').read_text()
        data=json.loads(re.search(r'<script id="recording" type="application/json">(.*?)</script>',text,re.S).group(1))
        data['fields']=unpack(data)
        return text,data
    def fixture(self,kind='dense'):
        capture=['phi','interface','deposit','activity','nutrient','chemical_potential'] if kind=='phase' else ['deposit','activity','nutrient']
        (self.root/'profile.toml').write_text(f'kind = "{kind}"\ntitle="Adapter fixture"\n[simulation]\nresolution=2\nseed=7\ncapture_fields='+json.dumps(capture)+'\n')
        run=dict(budget=dict(resolution=2,backend='metal',kind=kind,dt=0.123456789,extent=[2,2,2]),completed_steps=3)
        (self.root/'run.json').write_text(json.dumps(run))
        frames=[dict(prefix=f'frame-{step:08d}',step=step,scheduled_time=step*0.123456789,world_time=100+step,deposit=step*2) for step in [0,3]]
        (self.root/'frames.jsonl').write_text('\n'.join(json.dumps(f) for f in frames)+'\n')
        for f in frames:
            for field in capture:
                values=[-1+i/4+f['step']/10 for i in range(8)] if field=='phi' else [i/8+f['step']/10 for i in range(8)]
                self.write(f['prefix'],field,values)
        return frames
    def test_dense_scheduled_clock_metadata_and_quantization(self):
        frames=self.fixture();text,data=self.build()
        self.assertEqual(data['run']['shape'],[2,2,2]);self.assertEqual(data['run']['backend'],'metal')
        self.assertEqual([f['time'] for f in data['run']['frames']],[f['scheduled_time'] for f in frames])
        self.assertIn('METAL · 2³',text);self.assertIn('f.time.toFixed(6)',text)
        self.assertEqual(data['run']['config']['simulation']['dt'],0.123456789)
        lo,hi=data['scales']['deposit'];q=base64.b64decode(data['fields']['deposit'][1])
        for i,value in enumerate(q):self.assertLessEqual(abs(lo+value/255*(hi-lo)-(i/8+.3)),(hi-lo)/510+1e-6)
    def test_phase_uses_selected_signed_fields_not_auxiliary_files(self):
        frames=self.fixture('phase')
        # These are intentionally invalid: the viewer must not read unselected fields.
        for f in frames:
            for field in ['deposit','activity','chemical_potential']:(self.root/(f['prefix']+'-'+field+'.f32')).write_bytes(b'bad')
        text,data=self.build();self.assertEqual(data['run']['source_channels'],dict(deposit='phi',activity='interface',nutrient='nutrient'))
        self.assertLess(data['scales']['deposit'][0],0);self.assertIn('Slice / phi',text);self.assertIn('positive phase',text)
    def test_override_single_channel_marks_missing_channels(self):
        self.fixture('phase');text,data=self.build('--fields','phi')
        self.assertIsNone(data['run']['source_channels']['activity']);self.assertEqual(data['scales']['activity'],[0,0]);self.assertIn('disabled>Slice / not recorded',text)
        self.build('--fields','not_a_field',success=False)
    def test_wrong_byte_size_rejected_before_replacing_viewer(self):
        frames=self.fixture();(self.root/'viewer.html').write_text('previous viewer')
        (self.root/(frames[0]['prefix']+'-deposit.f32')).write_bytes(b'bad')
        result=self.build(success=False);self.assertIn('expected 32 bytes',result.stderr)
        self.assertEqual((self.root/'viewer.html').read_text(),'previous viewer')
    def test_legacy_constant_field(self):
        run=dict(shape=[2,2,2],backend='CPU',config={},frames=[dict(name='frame-0000',time=2.5,solid=0)])
        (self.root/'run.json').write_text(json.dumps(run))
        for field in ['deposit','activity','nutrient']:self.write('frame-0000',field,[.25]*8)
        text,data=self.build();self.assertEqual(data['scales']['deposit'],[.25,.25]);self.assertEqual(base64.b64decode(data['fields']['deposit'][0]),bytes(8))
        self.assertIn('max(ds.y-ds.x,1e-8)',text)
if __name__=='__main__':unittest.main()
