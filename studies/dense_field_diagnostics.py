"""Chunked native-grid diagnostics; never creates a resampled volume."""
import json,pathlib,sys
import numpy as np
root=pathlib.Path(sys.argv[1]);run=json.loads((root/'run.json').read_text());shape=run['shape'];count=int(np.prod(shape));cell_volume=float(np.prod(run['spacing']));report={'shape':shape,'extent':run['extent'],'cell_volume':cell_volume,'native_grid':True,'fields':{}}
for name in ('deposit','temperature','dose','pressure','intensity'):
    path=root/(name+'.f32')
    if not path.exists():continue
    assert path.stat().st_size==count*4
    field=np.memmap(path,dtype='<f4',mode='r',shape=(count,));total=0.;minimum=float('inf');maximum=-float('inf');finite=True;counts={str(v):0 for v in (.01,.1,.18,.5)}
    for start in range(0,count,1048576):
        chunk=field[start:start+1048576];total+=float(np.sum(chunk,dtype=np.float64));minimum=min(minimum,float(chunk.min()));maximum=max(maximum,float(chunk.max()));finite=finite and bool(np.isfinite(chunk).all())
        for threshold in counts:counts[threshold]+=int(np.count_nonzero(chunk>float(threshold)))
    report['fields'][name]={'finite':finite,'minimum':minimum,'maximum':maximum,'sum_float64':total,'physical_volume_integral':total*cell_volume,'fraction_above':{k:v/count for k,v in counts.items()}}
(root/'native-diagnostics.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
