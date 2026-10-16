#!/usr/bin/env python3
"""Compare recorded native fields and compact the exact energy ledger."""
import argparse,hashlib,json,pathlib
import numpy as np
p=argparse.ArgumentParser();p.add_argument('run',type=pathlib.Path);p.add_argument('dark',type=pathlib.Path);p.add_argument('--refined',type=pathlib.Path);a=p.parse_args();root=a.run;run=json.loads((root/'run.json').read_text());dark=json.loads((a.dark/'run.json').read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
assert run['shape']==dark['shape'];step=run['steps'];df=next(f for f in dark['frames'] if f['step']==step);cf=run['frames'][-1]
assert cf['time']==df['time'];assert run['source_sha256']['src/recurrent_material.jl']==dark['source_sha256']['src/recurrent_material.jl']
for name in ('solid','nutrient'):
    assert sha(root/run['frames'][0]['directory']/(name+'.f32'))==sha(a.dark/dark['frames'][0]['directory']/(name+'.f32'))
solid=np.fromfile(root/cf['directory']/'solid.f32',dtype='<f4').astype(float);control=np.fromfile(a.dark/df['directory']/'solid.f32',dtype='<f4').astype(float);dose=np.fromfile(root/cf['directory']/'dose.f32',dtype='<f4').astype(float);last=run['ledger'][-1]['power'];vv=(64/run['shape'][0])**3
summary=dict(native_shape=run['shape'],time=run['time'],steps=step,dt=run['dt'],elapsed_seconds=run['elapsed_seconds'],refresh_count=len(run['ledger']),rays_per_refresh=last['rays'],ray_ds=run['params']['ray_ds'],energy=run['energy'],relative_material_mass_residual=run['history'][-1]['mass_error']/run['initial_material_mass'],final_vs_dark_relative_solid=float(np.linalg.norm(solid-control)/np.linalg.norm(control)),same_initial_fields=True,same_recurrent_process_source=True,dark_run_sha256=sha(a.dark/'run.json'),raw_run_sha256=sha(root/'run.json'),final_dose_integral=float(dose.sum()*vv),final_deposited_power=last['deposited'],float32_dose_power_error=float(dose.sum()*vv-last['deposited']),maximum_unresolved_power=max(l['power']['unresolved'] for l in run['ledger']),power=[{k:v for k,v in l.items() if k!='endpoints'} for l in run['ledger']])
if a.refined:
    rr=json.loads((a.refined/'run.json').read_text());assert rr['shape']==run['shape'] and rr['time']==run['time']
    assert rr['source_sha256']['src/recurrent_material.jl']==run['source_sha256']['src/recurrent_material.jl']
    refsolid=np.fromfile(a.refined/rr['frames'][-1]['directory']/'solid.f32',dtype='<f4').astype(float)
    summary['native_half_ds_control']=dict(ray_ds=rr['params']['ray_ds'],raw_run_sha256=sha(a.refined/'run.json'),elapsed_seconds=rr['elapsed_seconds'],relative_final_solid=float(np.linalg.norm(refsolid-solid)/np.linalg.norm(solid)),relative_absorbed_energy=abs(rr['energy']['absorbed']-run['energy']['absorbed'])/run['energy']['absorbed'],energy=rr['energy'])
for key,path in [('pilot_mechanism_checks',pathlib.Path('outputs/art/illuminated-memory-gates/checks.json')),('native_ray_refinement',root/'native-ray-refinement.json')]:
    if path.exists():
        value=json.loads(path.read_text());summary[key]=value['metrics'] if key=='pilot_mechanism_checks' else value
(root/'feedback-summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k!='power'},indent=2))
