#!/usr/bin/env python3
"""Preserve measured lag traces and derive paired differences without state dumps."""
import hashlib,json,pathlib,sys

def compact(path):
    data=pathlib.Path(path).read_bytes();r=json.loads(data)
    baseline={(x['seed'],x['start']):x for x in r['results'] if x['amplitude']==0}
    out={k:v for k,v in r.items() if k!='results'}
    out['source_report_sha256']=hashlib.sha256(data).hexdigest()
    out['results']=[]
    max_residual=0.
    for x in r['results']:
        b=baseline[x['seed'],x['start']];paired=[]
        assert len(x['trace'])==len(b['trace'])
        for a,z in zip(x['trace'],b['trace']):
            assert a['offset']==z['offset']
            max_residual=max(max_residual,abs(a['material']),abs(a['potential']))
            delta={}
            for k in ('usable','uptake_a','uptake_b','oxidized','secretion_a','secretion_b','energy','enzyme','solubilized','work','heat','supplied','supplied_potential'):
                av,zv=a[k],z[k]
                delta[k]=[aa-zz if aa is not None and zz is not None else None for aa,zz in zip(av,zv)] if isinstance(av,list) else av-zv
            assert abs(delta['supplied'])<1e-10
            assert abs(delta['supplied_potential'])<1e-10
            paired.append({'offset':a['offset'],'delta':delta})
        out['results'].append({**x,'paired':paired})
    ratios=[x['feature_norm_squared'][i+4]/x['feature_norm_squared'][i] for x in r['contexts'] for i in range(4)]
    out['checks']={'max_abs_conservation_residual':max_residual,'feature_squared_norm_ratio_range':[min(ratios),max(ratios)],
        'all_residents_valid':all(all(t['resident_valid']) for x in r['results'] for t in x['trace']),
        'external_supply_matched':True}
    assert max_residual<1e-8
    return out
if __name__=='__main__':
    pathlib.Path(sys.argv[2]).write_text(json.dumps(compact(sys.argv[1]),separators=(',',':'),allow_nan=False)+'\n')
