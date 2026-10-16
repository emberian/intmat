#!/usr/bin/env python3
"""Compact actual matched allocation pulses without private controller tensors."""
import hashlib,json,statistics,sys
from pathlib import Path
source,out=map(Path,sys.argv[1:]);r=json.loads(source.read_text())
def observation(s):
    cells=s['cells']
    return {k:v for k,v in s.items() if k!='cells'}|{'population':len(cells),'resident_ids':[c['id'] for c in cells],
        'mean':{k:statistics.mean(c[k] for c in cells) for k in ('balance','total','trend','mass','energy','usable','allocation') if all(k in c for c in cells)}}
rows=[]
for row in r['rows']:
    rows.append({**row,'initial':observation(row['initial']),'samples':[
        {**s,'observation':observation(s['observation'])} for s in row['samples']]})
result={**r,'rows':rows,'raw_report_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
    'scope':'same-checkpoint matched physical branches; exact resident returns retained, controller tensors omitted'}
out.write_text(json.dumps(result,indent=2)+'\n')
