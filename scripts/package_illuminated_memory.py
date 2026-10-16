#!/usr/bin/env python3
"""Prepare compact publication inputs without altering the raw coupled epoch."""
import argparse,hashlib,json,pathlib,shutil
p=argparse.ArgumentParser();p.add_argument('run',type=pathlib.Path);a=p.parse_args();root=a.run;out=root/'publication';out.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
run=json.loads((root/'run.json').read_text());ledger=run.pop('ledger');run['raw_run_sha256']=sha(root/'run.json');run['illumination_summary']=[{k:v for k,v in item.items() if k!='endpoints'} for item in ledger];run['raw_endpoint_ledger_retained']=True
(out/'run.json').write_text(json.dumps(run,indent=2))
for name in ('final.png','final.png.optical.json','render-job.json','recipe.toml','final-status.json','art-reviewed.json','viewer.html'):
    shutil.copyfile(root/name,out/name)
for name in ('feedback-summary.json','native-ray-refinement.json'):
    if (root/name).exists():shutil.copyfile(root/name,out/name)
print(out)
