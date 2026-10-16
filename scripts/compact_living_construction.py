#!/usr/bin/env python3
"""Retain measured living-world identities, budgets and events without raw fields."""
import hashlib,json,math,pathlib,sys

def compact(path):
    data=pathlib.Path(path).read_bytes();r=json.loads(data);frames=[]
    for f in r['frames']:
        assert max(abs(f['ledger'][k]) for k in ('material_residual','potential_residual'))<1e-8
        frames.append({**{k:v for k,v in f.items() if k!='material'},'material_totals':{k:sum(v) for k,v in f['material'].items()}})
    founders={c['id']:c for c in frames[0]['cells']};extent=r['metadata']['extent'];motion=[]
    for c in frames[-1]['cells']:
        if c['id'] in founders:
            motion.append(math.sqrt(sum(min(abs(a-b),e-abs(a-b))**2 for a,b,e in zip(c['position'],founders[c['id']]['position'],extent))))
    return {**{k:v for k,v in r.items() if k!='frames'},'frames':frames,'source_report_sha256':hashlib.sha256(data).hexdigest(),
        'summary':{'surviving_founder_mean_displacement':sum(motion)/len(motion) if motion else None,
          'max_generation':max(c['generation'] for f in frames for c in f['cells']),
          'paid_nonzero_copies':sum(e['copy_funded'] and e['copy_work']>0 for e in r['birth_events']),
          'copy_work':sum(e['copy_work'] for e in r['birth_events']),
          'max_abs_conservation_residual':max(abs(f['ledger'][k]) for f in frames for k in ('material_residual','potential_residual'))}}
if __name__=='__main__':
    pathlib.Path(sys.argv[1]).write_text(json.dumps({'schema':'living-construction-comparison-v1','runs':[compact(p) for p in sys.argv[2:]]},separators=(',',':'),allow_nan=False)+'\n')
