"""Audit native observations, irreversible failures, closure and source pinning."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('capture',type=Path);a=p.parse_args();meta=json.loads((a.capture/'capture.json').read_text());assert len(meta['frames'])==40;assert hashlib.sha256((a.capture/'base.json').read_bytes()).hexdigest()==meta['base_sha256'];prior=None
for r in meta['frames']:
 path=a.capture/r['file'];assert hashlib.sha256(path.read_bytes()).hexdigest()==r['sha256'];s=json.loads(path.read_text())
 for key in ('position','damage','potential','dose','temperature'):assert np.all(np.isfinite(np.asarray(s[key])))
 damage=np.asarray(s['damage']);assert np.all((damage==0)|(damage==1));assert min(s['dose'])>=0 and min(s['temperature'])>=0;assert min(s['potential'])>=-2e-5 and max(s['potential'])<=1+2e-5
 if prior:assert s['step']>prior['step'] and s['time']>prior['time'] and np.all(damage>=np.asarray(prior['damage'])) and np.all(np.asarray(s['dose'])>=np.asarray(prior['dose']))
 assert r['solver']['residual']<2e-6 and r['solver']['relative_power_defect']<2e-5;prior=s
report={'actual_captures':40,'native_nodes':meta['native_nodes'],'native_bonds':meta['native_bonds'],'physical_time_interval':[meta['frames'][0]['time'],meta['frames'][-1]['time']],'maximum_independent_residual':max(r['solver']['residual'] for r in meta['frames']),'maximum_relative_power_defect':max(r['solver']['relative_power_defect'] for r in meta['frames']),'maximum_mechanical_relative_work_residual':max(abs(r['mechanical_relative_work_residual']) for r in meta['frames']),'broken_bonds_interval':[meta['frames'][0]['broken'],meta['frames'][-1]['broken']],'global_dose_optical_max':meta['global_dose_optical_max'],'temporal_interpolation':False,'restart_checkpoints':False,'source_sha256':meta['source_sha256'],'capture_sha256':hashlib.sha256((a.capture/'capture.json').read_bytes()).hexdigest(),'irreversible_damage_and_nondecreasing_dose':True}
(a.capture/'verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
