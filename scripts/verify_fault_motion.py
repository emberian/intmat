"""Check the fresh stored trajectory and its exact approved-final-state boundary."""
import argparse,hashlib,json,math
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('capture',type=Path);p.add_argument('approved',type=Path);a=p.parse_args();files=sorted((a.capture/'states').glob('*.json'));assert len(files)==40
baseline=json.loads(a.approved.read_text());prior=None;records=[]
for path in files:
 raw=path.read_bytes();s=json.loads(raw)
 assert all(math.isfinite(x) for key in ['position','reference','potential','temperature','dose','strain'] for x in s[key])
 assert min(s['dose'])>=0 and min(s['temperature'])>=0 and min(s['potential'])>=0 and max(s['potential'])<=1
 if prior:
  assert s['step']>prior['step']
  assert all(a>=b for a,b in zip(s['damage'],prior['damage']))
  assert all(a>=b for a,b in zip(s['dose'],prior['dose']))
 records.append({'step':s['step'],'time':s['time'],'sha256':hashlib.sha256(raw).hexdigest(),'broken':sum(s['damage']),'maximum_dose':max(s['dose'])});prior=s
keys=['dims','position','reference','bonds','damage','strain','dose','temperature','potential']
assert all(prior[k]==baseline[k] for k in keys), 'Captured final state differs from approved still'
report={'frames':records,'actual_stored_states':True,'interpolation':False,'final_state_exact_match':keys,'approved_state_sha256':hashlib.sha256(a.approved.read_bytes()).hexdigest(),'finite_fields':True,'irreversible_damage':True,'nondecreasing_dose':True,'potential_in_unit_interval':True,'source_run':json.loads((a.capture/'run.json').read_text())}
(a.capture/'capture-verification.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k not in ['frames','source_run']},indent=2))
