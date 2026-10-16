#!/usr/bin/env python3
"""Keep energetic and developmental returns distinct, with fixed denominators."""
import hashlib,json,statistics,sys
from pathlib import Path
source,out=map(Path,sys.argv[1:]);r=json.loads(source.read_text());arms={}
for name,a in r['arms'].items():
    probes={}
    for p in a['probes']:
        values=p['returns']
        probes[p['mode']]={'initial_count':len(values),'survivors':sum(x['alive'] for x in values),
            'mean_inventory_return':statistics.mean(x['all_initial_value'] for x in values),
            'mean_development_return':statistics.mean(x['development_value'] for x in values),
            'partition_mass':sum(x['partition_mass'] for x in values)}
    samples=[{k:v for k,v in s.items() if k!='cells'}|{'population':len(s['cells']),
        'mean_usable':statistics.mean(c['usable'] for c in s['cells']),
        'mean_development_value':statistics.mean(c['development_value'] for c in s['cells']),
        'biomass':sum(c['mass'] for c in s['cells'])} for s in a['samples']]
    arms[name]={k:v for k,v in a.items() if k not in ('samples','probes','events','capital_gifts')}|{'samples':samples,'probes':probes,'birth_count':len(a['events'])}
result={k:v for k,v in r.items() if k!='arms'}|{'arms':arms,'raw_report_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
    'scope':'same-arm physical-state frozen interventions; dead endpoint0 in both values; actual same-parent child U/mass gifts; common noise is not independent replication'}
out.write_text(json.dumps(result,indent=2)+'\n')
