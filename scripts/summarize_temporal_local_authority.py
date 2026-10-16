#!/usr/bin/env python3
import hashlib,json,math,statistics,sys
from pathlib import Path
source,out=map(Path,sys.argv[1:]);r=json.loads(source.read_text());neutral=r['results'][0];rows=[]
for branch in r['results'][1:]:
    samples=[]
    for s,n in zip(branch['samples'],neutral['samples']):
        base={x['id']:x['all_initial_value'] for x in n['returns']}
        diff={x['id']:x['all_initial_value']-base[x['id']] for x in s['returns']}
        others=[v for k,v in diff.items() if k!=branch['producer']]
        cells={c['id']:c for c in s['observation']['cells']};ncells={c['id']:c for c in n['observation']['cells']}
        own=None
        if branch['producer'] in cells and branch['producer'] in ncells:
            a=cells[branch['producer']];b=ncells[branch['producer']]
            own={k:a[k]-b[k] for k in ('mass','energy','uptake_a','uptake_b')}
            own.update(return_difference=diff[branch['producer']],neutral_return=base[branch['producer']])
        globals={k:[a-b for a,b in zip(s['observation'][k],n['observation'][k])] for k in ('catalyst','bound','converted','funded','requested')}
        samples.append({'lag':s['lag'],'own':own,'other_count':len(others),'others_sum_return_difference':sum(others),
            'others_rms_return_difference':math.sqrt(statistics.mean(v*v for v in others)),
            'initial_count':len(s['returns']),'survivors':sum(x['alive'] for x in s['returns']),
            'neutral_survivors':sum(x['alive'] for x in n['returns']),'global_differences':globals,
            'material_residual':s['observation']['ledger']['material_residual'],'potential_residual':s['observation']['ledger']['potential_residual']})
    rows.append({k:v for k,v in branch.items() if k!='samples'}|{'samples':samples})
result={k:v for k,v in r.items() if k not in ('results','initial')}|{'results':rows,'raw_report_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
out.write_text(json.dumps(result,indent=2)+'\n')
