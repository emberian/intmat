#!/usr/bin/env python3
import array,hashlib,json,sys,tempfile,unittest
from pathlib import Path
from summarize_native_observatory import summarize
class GraphTests(unittest.TestCase):
    def test_reciprocal_only_and_pinned_energy_secretion(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);binary=bytearray();specs={}
            for name,kind,values,shape in [('ids','i',[1,2,3,4],[4]),('bond_ids','i',[2,1,2,99,0,3,4,0],[4,2]),('energy','d',[1,2,3,4],[4]),('last_step_secretion','f',[0,1,0,0,0,0,2,0],[4,2])]:
                data=array.array(kind,values).tobytes();specs[name]={'offset':len(binary),'bytes':len(data),'type':{'i':'i32','d':'f64','f':'f32'}[kind],'shape':shape};binary.extend(data)
            (root/'state.bin').write_bytes(binary);m=dict(arrays=specs,binary='state.bin',binary_sha256=hashlib.sha256(binary).hexdigest(),endian=sys.byteorder,time=0,step=0,shape=[8,8,8],source_sha256='source')
            (root/'frame.json').write_text(json.dumps(m));(root/'observatory.json').write_text(json.dumps({'frames':[{'file':'frame.json','time':0,'side':8}]}))
            f=summarize(root)['frames'][0];self.assertEqual(f['component_sizes'],[3,1]);self.assertEqual(f['reciprocal_edges'],2);self.assertEqual(f['nonreciprocal_records'],[[3,4]]);self.assertEqual(f['invalid_records'],[[4,99]]);self.assertEqual(f['mean_energy'],2.5);self.assertEqual(f['last_step_secretors'],2);self.assertEqual(f['last_step_secretion'],[1,2])
            (root/'state.bin').write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError,'SHA'):summarize(root)
if __name__=='__main__':unittest.main()
