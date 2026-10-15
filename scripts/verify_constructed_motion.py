#!/usr/bin/env python3
import argparse,array,base64,gzip,json,math,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description='Verify the compact exhibit against every original observation and exact paid-flow sidecar.')
parser.add_argument('--source',type=Path,default=root/'outputs/constructed-feedback-observation-20261015')
parser.add_argument('--exhibit',type=Path,default=root/'outputs/constructed-motion-20261015')
parser.add_argument('--report',type=Path)
args=parser.parse_args();source=args.source;out=args.exhibit
data=json.loads(gzip.decompress((out/'recording.json.gz').read_bytes())); counts={'frames':0,'cells':0,'field_samples':0,'bond_endpoints':0,'paid_values':0}; maximum={'geometry_error_physical':0.,'neural_error':0.,'paid_relative_error':0.}
def unpack(kind,text):
 a=array.array(kind);a.frombytes(base64.b64decode(text))
 if sys.byteorder!='little':a.byteswap()
 return list(a)
for run in data['runs']:
 previous={key:[0]*1024 for key in data['fields']}; flows={r['step']:r for r in map(json.loads,(source/'analysis'/f"paid-flow-{run['seed']}"/'interval-cell-flows.jsonl').read_text().splitlines())}
 for frame in run['frames']:
  original=json.loads((source/f"strong-{run['seed']}-continuation-3100"/f"frame-{round(frame['time']*10):07d}.json").read_text()); ids=unpack('H',frame['ids']); geometry=unpack('H',frame['geometry']); neural=unpack('b',frame['neural']); paid=unpack('f',frame['paid']); flow={c['id']:c for c in flows[original['step']]['cells']}
  assert [run['identities'][i][0] for i in ids]==[c['id'] for c in original['agents']]
  for key in data['fields']:
   restored=[(a+b)%256 for a,b in zip(unpack('B',frame['fields'][key]),previous[key])]; expected=[round(v/(data['scales'][key] or 1)*255) for v in original['fields'][key]['values']]; assert restored==expected;previous[key]=restored;counts['field_samples']+=len(restored)
  bonds=unpack('H',frame['bonds']);assert [run['identities'][i][0] for i in bonds]==[i for edge in original['metrics']['cell_bonds'] for i in edge];counts['bond_endpoints']+=len(bonds)
  for i,c in enumerate(original['agents']):
   for d in range(2):
    restored=geometry[i*3+d]/65535*run['grid']['extent'][d]; expected=c['cell_position'][d]/run['grid']['shape'][d]*run['grid']['extent'][d];error=abs(restored-expected);assert error<=run['grid']['extent'][d]/65535/2+1e-12; maximum['geometry_error_physical']=max(maximum['geometry_error_physical'],error)
   original_neural=c['hidden']+c['adaptation']+c['actions']
   for a,b in zip(neural[i*27:i*27+27],original_neural):
    error=abs(a/127-b);assert error<=1/254+1e-9;maximum['neural_error']=max(maximum['neural_error'],error)
   row=flow.get(c['id'])
   if row:
    expected=row['actual_secretion']+[row['actual_construction_work'],row['actual_motor_work']]
    for a,b in zip(paid[i*4:i*4+4],expected):
     assert abs(a-b)<=abs(b)*1e-7+1e-44
     if b:maximum['paid_relative_error']=max(maximum['paid_relative_error'],abs(a-b)/abs(b))
     counts['paid_values']+=1
   else:assert all(math.isnan(v) for v in paid[i*4:i*4+4])
  counts['frames']+=1;counts['cells']+=len(ids)
report={'passed':True,'checks':counts,'maximum_observed_encoding_errors':maximum,'interpretation':'All displayed IDs/bonds match recorded state; all field bytes equal direct source quantization; all geometry/neural/paid-flow errors stay within declared encoding bounds.'}
if args.report:
 args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
