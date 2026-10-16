#!/usr/bin/env python3
"""Compact measured continuous-policy outcomes; preserve survival denominators."""
import hashlib,json,statistics,sys
from pathlib import Path
source,out=map(Path,sys.argv[1:])
r=json.loads(source.read_text());arms={}
for name,a in r['arms'].items():
    probes={}
    for p in a['probes']:
        initial={c['id']:c for c in p['initial']['cells']}
        alive=[x['value'] for x in p['returns'] if x['alive']]
        # New reports retain exact parent gifts even when the parent dies.
        # Legacy records can safely supply a whole-cohort mean only without deaths.
        probes[p['mode']]={'initial_count':len(initial),'survivors':len(alive),
            'survivor_conditional_mean_return':statistics.mean(alive) if alive else None,
            'all_initial_mean_return':statistics.mean(x['all_initial_value'] for x in p['returns']) if all('all_initial_value' in x for x in p['returns']) else (statistics.mean(x['value'] for x in p['returns']) if len(alive)==len(initial) else None),
            'all_initial_metric_note':('dead endpoint usable=0, plus actual same-parent partition gifts' if all('all_initial_value' in x for x in p['returns']) else 'legacy report: whole-cohort mean withheld if any death')}
    samples=[]
    for s in a['samples']:
        samples.append({k:v for k,v in s.items() if k!='cells'}|{'population':len(s['cells']),
            'mean_usable':statistics.mean(c['usable'] for c in s['cells']) if s['cells'] else None})
    arms[name]={k:v for k,v in a.items() if k not in ('samples','probes','events')}|{'samples':samples,'probes':probes,'birth_count':len(a['events'])}
result={k:v for k,v in r.items() if k!='arms'}|{'arms':arms,'raw_report_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
    'summary_scope':'same-arm physical-state frozen comparisons; matched noise across arms; no independent-noise statistical claim'}
out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n')
