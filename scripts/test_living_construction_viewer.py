"""Schema, identity and x-fast field checks for the actual-capture adapter."""
import base64
import copy
import gzip
import json
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from build_living_construction_viewer import FIELDS, LEDGERS, build, pack, validate_packed

def fixture():
    def cell(identity,parent=0,generation=0):
        return dict(id=identity,parent_id=parent,lineage=1,generation=generation,position=[.5+.3*identity,1,1],radius=.11,mass=.2,energy=.3,stores=[.1,.1],learned_l2=.002,learning_work=.0001)
    cells=[[cell(1)],[cell(1),cell(2,1,1)],[cell(2,1,1)]]
    frames=[]
    for j,time in enumerate([0,.5,2]):
        fields={key:[100*j+10*k+i/100 for i in range(24)] for k,key in enumerate(FIELDS)}
        frames.append(dict(time=time,cells=cells[j],material=fields,ledger={k:0 for k in LEDGERS},births=int(j>0),deaths=int(j==2)))
    return dict(schema='living-construction-v1',synthetic=True,metadata=dict(shape=[4,3,2],extent=[2,2,2],origin=[0,0,0],dt=.02,seed=211,source_rates=[.4,.4],source_budgets=[12.8,12.8],horizon=2,actor_rate=.003,sigma=.3,inheritance_fraction=.5,birth_transfer_adjusted=True,material_units='amount per voxel',field_order='x-fast'),frames=frames,birth_events=[dict(time=.5,parent_id=1,child_id=2,child_usable_before_copy=.6,copied_l1=.02,copy_work=.0002,copy_funded=True,adjusted_pending_count=2)])

class LivingViewerTest(unittest.TestCase):
    def test_actual_time_identity_and_field_mapping(self):
        p=pack(fixture(),'a'*64);raw=gzip.decompress(base64.b64decode(p['fields']['data']));values=struct.unpack('<'+str(len(raw)//4)+'f',raw)
        self.assertEqual([f['time'] for f in p['frames']],[0,.5,2])
        self.assertEqual([c['id'] for c in p['frames'][2]['cells']],[2])
        self.assertEqual(p['frames'][1]['cells'][1]['parent_id'],1)
        self.assertEqual(p['frames'][2]['field_offset'],2*8*24)
        # x=3,y=1,z=1: x + nx*(y + ny*z) = 19, channel=5, frame=2.
        self.assertAlmostEqual(values[2*8*24+5*24+19],250.19,places=4)
        self.assertEqual(p['fields']['maximum'][5],250.23)
        self.assertTrue(p['synthetic']);self.assertEqual(p['birth_events'][0]['adjusted_pending_count'],2)
    def test_duplicate_id_rejected(self):
        f=fixture();f['frames'][0]['cells']*=2
        with self.assertRaisesRegex(ValueError,'Duplicate'):pack(f,'a'*64)
    def test_reused_identity_rejected(self):
        f=fixture();f['frames'][1]['cells'][0]['lineage']=2
        with self.assertRaisesRegex(ValueError,'ancestry'):pack(f,'a'*64)
    def test_invalid_field_and_time_rejected(self):
        for mode in ('field','time','order'):
            f=fixture()
            if mode=='field':f['frames'][0]['material']['raw_a'][0]=float('nan')
            elif mode=='time':f['frames'][2]['time']=-1
            else:f['metadata']['field_order']='z-fast'
            with self.assertRaises(ValueError):pack(f,'a'*64)
    def test_bad_radius_and_event_rejected(self):
        f=fixture();f['frames'][0]['cells'][0]['radius']=-1
        with self.assertRaisesRegex(ValueError,'radius'):pack(f,'a'*64)
        f=fixture();f['birth_events'][0]['copy_funded']='yes'
        with self.assertRaisesRegex(ValueError,'Boolean'):pack(f,'a'*64)
    def test_packed_rebuild_and_corruption(self):
        p=pack(fixture(),'a'*64);self.assertEqual(validate_packed(p),p)
        for change in ('hash','offset','bytes','scale'):
            q=copy.deepcopy(p)
            if change=='hash':q['fields']['sha256']='0'*64
            elif change=='offset':q['frames'][1]['field_offset']=0
            elif change=='bytes':q['fields']['bytes']+=4
            else:q['fields']['maximum'][0]+=1
            with self.assertRaises(ValueError):validate_packed(q)
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);source=root/'raw.json';source.write_text(json.dumps(fixture()))
            build(source,root/'first',packed_output=root/'packed.json');build(root/'packed.json',root/'second')
            self.assertEqual((root/'first/index.html').read_bytes(),(root/'second/index.html').read_bytes())
    def test_no_silent_artifact_budget_expansion(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);source=root/'source.json';source.write_text(json.dumps(fixture()))
            with self.assertRaisesRegex(ValueError,'above explicit'):build(source,root/'viewer',100)
            self.assertFalse((root/'viewer/index.html').exists())
            p,size=build(source,root/'viewer');self.assertLess(size,1_000_000)
            self.assertIn('Synthetic private fixture',(root/'viewer/index.html').read_text())

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='--fixture':
        root=Path(sys.argv[2]);root.mkdir(parents=True,exist_ok=True);source=root/'synthetic.json';source.write_text(json.dumps(fixture()));build(source,root)
    else:unittest.main()
