#!/usr/bin/env python3
"""Read-only SHA-pinned graph/energy/final-step secretion summaries."""
import argparse,array,hashlib,json,sys
from pathlib import Path

def summarize(root):
    root=Path(root);raw=(root/'observatory.json').read_bytes();index=json.loads(raw);rows=[];previous=None
    best={}
    for entry in index['frames']:
        if entry['time'] not in best or entry['side']>best[entry['time']]['side']:best[entry['time']]=entry
    for entry in sorted(best.values(),key=lambda f:f['time']):
        path=root/entry['file'];metadata_raw=path.read_bytes();m=json.loads(metadata_raw);binary=(path.parent/m['binary']).read_bytes()
        if hashlib.sha256(binary).hexdigest()!=m['binary_sha256']:raise ValueError('binary SHA mismatch')
        def values(name):
            s=m['arrays'][name];v=array.array({'i32':'i','f32':'f','f64':'d'}[s['type']]);v.frombytes(binary[s['offset']:s['offset']+s['bytes']])
            if m['endian']!=sys.byteorder:v.byteswap()
            return v
        ids=list(values('ids'));n=len(ids);alive=set(ids);slots=values('bond_ids');directed={(ids[i],slots[i+n*k]) for i in range(n) for k in range(m['arrays']['bond_ids']['shape'][1]) if slots[i+n*k]>0}
        foreign=sorted([list(e) for e in directed if e[1] not in alive or e[0]==e[1]])
        unmatched=sorted([list(e) for e in directed if e[1] in alive and e[0]!=e[1] and (e[1],e[0]) not in directed])
        edges={(a,b) for a,b in directed if a<b and (b,a) in directed and b in alive};adj={i:set() for i in ids}
        for a,b in edges:adj[a].add(b);adj[b].add(a)
        components=[];remaining=set(ids)
        while remaining:
            stack=[min(remaining)];remaining.remove(stack[0]);component=[]
            while stack:
                i=stack.pop();component.append(i)
                for j in adj[i]&remaining:remaining.remove(j);stack.append(j)
            components.append(sorted(component))
        energy=values('energy');secretion=values('last_step_secretion');union=edges|previous if previous is not None else None
        rows.append(dict(time=m['time'],step=m['step'],side=m['shape'][0],population=n,source_sha256=m['source_sha256'],metadata_sha256=hashlib.sha256(metadata_raw).hexdigest(),binary_sha256=m['binary_sha256'],reciprocal_edges=len(edges),component_sizes=sorted(map(len,components),reverse=True),components=components,nonreciprocal_records=unmatched,invalid_records=foreign,edge_jaccard_previous=None if union is None else len(edges&previous)/len(union) if union else 1,mean_energy=sum(energy)/n if n else None,last_step_secretors=sum(secretion[i]>0 or secretion[i+n]>0 for i in range(n)),last_step_secretion=[sum(secretion[:n]),sum(secretion[n:])]))
        previous=edges
    return dict(schema='native-observatory-graph-summary-1',manifest_sha256=hashlib.sha256(raw).hexdigest(),note='Reciprocal live-ID bonds only; components describe bond-connected cells. Secretion amounts cover the final physical step, not the capture interval.',frames=rows)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('observatory',type=Path);p.add_argument('output',type=Path);a=p.parse_args();result=summarize(a.observatory)
    # Preserve earlier summaries; caller chooses a new evidence path each time.
    with a.output.open('x') as f:json.dump(result,f,separators=(',',':'));f.write('\n')
    print(json.dumps({'frames':len(result['frames']),'output':str(a.output)}))
if __name__=='__main__':main()
