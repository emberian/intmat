#!/usr/bin/env python3
"""Publish observed joint-policy traces, excluding learned tensors and local paths."""
import hashlib,json,math,pathlib,sys

def mean(values):
    return sum(values)/len(values) if values else None

def compact(path):
    path=pathlib.Path(path);data=path.read_bytes();r=json.loads(data)
    out={k:v for k,v in r.items() if k!='results'}
    out['source_report_sha256']=hashlib.sha256(data).hexdigest();out['results']={}
    for arm,x in r['results'].items():
        summaries={}
        for context in ('continued','shifted'):
            learned=x['branches']['learned_'+context];erased=x['branches']['erased_'+context]
            deltas=[a-b if a is not None and b is not None else None for a,b in zip(learned['return_value'],erased['return_value'])]
            if arm=='no_credit':
                assert all(v==0 for v in deltas if v is not None)
            expr=learned['expression'];reserve=learned['reserve_fraction']
            def logit(v): return math.log(v/(1-v)) if 0<v<1 else None
            offsets={
                'investment':[logit(v)+.8 if logit(v) is not None else None for v in expr[:4]],
                'allocation':[logit(v) for v in expr[4:8]],
                'reserve':[.5*logit(v)+math.atanh(.6) if logit(v) is not None else None for v in reserve]}
            summaries[context]={'cell_delta':deltas,'mean_delta':mean([v for v in deltas if v is not None]),
                'actual_final_held_logit_displacement':offsets}
        out['results'][arm]={'trace':x['trace'],'branches':x['branches'],'frozen_comparison':summaries}
    return out

if __name__=='__main__':
    output=pathlib.Path(sys.argv[1]);reports=[compact(p) for p in sys.argv[2:]]
    output.write_text(json.dumps({'schema':'native-joint-construction-comparison-v1','runs':reports},separators=(',',':'),allow_nan=False)+'\n')
