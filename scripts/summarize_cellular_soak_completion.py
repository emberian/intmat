#!/usr/bin/env python3
"""Read-only compact longitudinal summary of completed cellular soaks."""
import argparse
import hashlib
import json
from pathlib import Path
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize(root):
    root = Path(root)
    rows = []
    chunks = []
    for path in sorted(root.glob('metrics-*.jsonl')):
        raw = path.read_bytes()
        if raw and not raw.endswith(b'\n'):
            raise ValueError(f'Incomplete metrics tail: {path}')
        values = [json.loads(line) for line in raw.splitlines()]
        rows.extend(values)
        chunks.append(dict(file=path.name, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest(), rows=len(values)))
    if not rows:
        raise ValueError(f'No metrics: {root}')
    if any(a['step'] >= b['step'] for a,b in zip(rows,rows[1:])):
        raise ValueError('Expected strictly increasing steps in one completed session')
    keys = ('population', 'births', 'deaths', 'bonds', 'max_generation', 'surviving_lineages',
            'biomass', 'mean_age', 'junction_transferred', 'material_error', 'chemical_energy_error',
            'blocked_births', 'synthesized', 'oxidized', 'chemical_heat', 'work_spent')
    def compact(row):
        return {**{k: row[k] for k in ('step','time','observed_at_unix')},
                'metrics': {k:row['metrics'][k] for k in keys if k in row['metrics']}}
    extrema = {}
    for key in keys:
        values = [r for r in rows if isinstance(r['metrics'].get(key), (int,float))]
        if values:
            low,high = min(values,key=lambda r:r['metrics'][key]),max(values,key=lambda r:r['metrics'][key])
            extrema[key] = {'min':low['metrics'][key], 'min_step':low['step'],
                            'max':high['metrics'][key], 'max_step':high['step']}
    # One sampled point per ten wall-clock minutes, plus the exact first/final row.
    selected = [rows[0]]
    for row in rows[1:-1]:
        if row['observed_at_unix']-selected[-1]['observed_at_unix'] >= 600:
            selected.append(row)
    if selected[-1] is not rows[-1]:
        selected.append(rows[-1])
    bonded = [r for r in rows if r['metrics'].get('bonds',0)>0]
    lineage_index = next((i for i,r in enumerate(rows) if r['metrics'].get('surviving_lineages')==1),None)
    pointer=json.loads((root/'checkpoint.json').read_text())
    cp=root/pointer['file']
    before=cp.stat();sha=digest(cp);after=cp.stat()
    if (before.st_ino,before.st_size,before.st_mtime_ns)!=(after.st_ino,after.st_size,after.st_mtime_ns):
        raise ValueError('Checkpoint changed during summary')
    archives=sorted(root.glob('state-*.jls'))
    status=json.loads((root/'status.json').read_text())
    if status['state']!='stopped' or status['reason']!='wall_time_complete':
        raise ValueError(f'Not a normally completed soak: {root}')
    source=root.parent/'source';provenance=json.loads((source/'provenance.json').read_text())
    mismatch=[name for name,sha in provenance['hashes'].items() if digest(source/name)!=sha]
    return {'run':str(root),'status':status,'metrics_rows':len(rows),'metrics_chunks':chunks,
            'observed_wall_seconds':rows[-1]['observed_at_unix']-rows[0]['observed_at_unix'],
            'first':compact(rows[0]),'final':compact(rows[-1]),'extrema':extrema,
            'sampled_trajectory':list(map(compact,selected)),
            'positive_bond_rows':sum(r['metrics'].get('bonds',0)>0 for r in rows),
            'last_observed_positive_bonds':compact(bonded[-1]) if bonded else None,
            'single_lineage_transition_bracket':[compact(rows[max(0,lineage_index-1)]),compact(rows[lineage_index])] if lineage_index is not None else None,
            'first_observed_single_lineage':next((compact(r) for r in rows if r['metrics'].get('surviving_lineages')==1),None),
            'first_observed_zero_population':next((compact(r) for r in rows if r['metrics']['population']==0),None),
            'checkpoint':{'pointer':pointer,'bytes':after.st_size,'sha256':sha,'stable':True},
            'archives':{'count':len(archives),'bytes':sum(p.stat().st_size for p in archives),
                        'first':archives[0].name if archives else None,'last':archives[-1].name if archives else None},
            'source_verification':{'files':len(provenance['hashes']),'mismatches':mismatch,
                                   'manifest_sha256':digest(source/'provenance.json')}}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs',nargs='+')
    args=parser.parse_args()
    print(json.dumps({'checked_unix':time.time(),'method':'All complete metrics rows summarized; trajectory sampled every >=600 wall seconds; no checkpoint deserialization or simulation steps',
                      'runs':[summarize(p) for p in args.runs]},indent=2))
