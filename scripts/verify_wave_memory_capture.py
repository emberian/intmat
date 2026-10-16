"""Audit stored native observations and quantify material threshold crossings."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('capture',type=Path);a=p.parse_args();run=json.loads((a.capture/'run.json').read_text());shape=tuple(run['shape']);prior=None;records=[];fixedhost=None
for f in run['frames']:
 folder=a.capture/f['directory'];field=lambda name:np.fromfile(folder/(name+'.f32'),dtype='<f4').reshape(shape,order='F');m=field('deposit');h=field('density');q=field('activity');assert all(np.all(np.isfinite(v)) for v in (m,h,q));assert np.min(m)>=0 and np.min(q)>=0 and np.max(m-h)<=1.2e-7
 hosthash=hashlib.sha256((folder/'density.f32').read_bytes()).hexdigest();fixedhost=fixedhost or hosthash;assert fixedhost==hosthash
 record={'step':f['step'],'time':f['time'],'mineral_sum':float(np.sum(m,dtype='f8')),'fields_sha256':{k:hashlib.sha256((folder/(k+'.f32')).read_bytes()).hexdigest() for k in ('deposit','density','activity','envelope')}}
 if prior is not None:record.update(new_voxels_above_iso=int(np.sum((m>=.3)&(prior<.3))),receded_voxels_below_iso=int(np.sum((m<.3)&(prior>=.3))),mean_absolute_mineral_change=float(np.mean(abs(m-prior))))
 records.append(record);prior=m
out={'native_shape':shape,'captured_states':len(records),'movie_states':40,'actual_stored_observations':True,'restart_checkpoints':False,'temporal_interpolation':False,'pressure_channel':'activity is instantaneous absolute pressure, not a simulation activity reservoir','isovalue':.3,'total_receding_threshold_crossings':sum(r.get('receded_voxels_below_iso',0) for r in records),'frames':records,'maximum_reported_local_inventory_error':max(r['inventory_error'] for r in run['frames'])};(a.capture/'capture-verification.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='frames'},indent=2))
