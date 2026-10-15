#!/usr/bin/env python3
"""Plot detached, actual GPU pilot fields and selected individual histories."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('export',type=Path);p.add_argument('output',type=Path)
a=p.parse_args();a.output.mkdir()
m=json.loads((a.export/'frames.json').read_text());frames=m['frames']
def field(f):
    n=f['coarse'];return np.fromfile(a.export/(f['name']+'.fields.f32'),dtype='<f4').reshape((n,n,n,9),order='F')
def array(f,name,dtype,shape):
    return np.fromfile(a.export/(f['name']+'.'+name),dtype=dtype).reshape(shape,order='F')
fields=[field(f) for f in frames];t=np.array([f['time'] for f in frames]);last=frames[-1];L=last['extent'][0]
fig,ax=plt.subplots(2,2,figsize=(11,8.2),constrained_layout=True)
cat=(fields[-1][...,4]+fields[-1][...,5]).sum(axis=2)*(L/last['coarse'])
im=ax[0,0].imshow(cat.T,origin='lower',extent=(0,L,0,L),cmap='magma');fig.colorbar(im,ax=ax[0,0],label='mass / area')
ax[0,0].set(title='Catalyst present: column integral',xlabel='x',ylabel='y')
delta=(fields[-1][...,8]-fields[0][...,8]).mean(axis=2);v=float(np.max(np.abs(delta)))
im=ax[0,1].imshow(delta.T,origin='lower',extent=(0,L,0,L),cmap='RdBu_r',vmin=-v,vmax=v);fig.colorbar(im,ax=ax[0,1],label='final − initial')
ax[0,1].set(title='Combined permeability: column-mean change',xlabel='x',ylabel='y')
for identity in (1,2,3,4):
    values=[]
    for f in frames:
        ids=array(f,'trace_ids.i32','<i4',(f['trace_count'],));h=array(f,'hidden.f32','<f4',(f['trace_count'],f['hidden']))
        loc=np.flatnonzero(ids==identity);values.append(float(h[loc[0],0]) if len(loc) else np.nan)
    ax[1,0].plot(t,values,'o-',ms=4,label=f'ID {identity}')
ax[1,0].set(title='Neural channel 1 in four actual cells',xlabel='time',ylabel='hidden state');ax[1,0].legend(frameon=False,ncol=2)
volume=(L/last['coarse'])**3
catalyst=np.array([np.sum(f[...,4:6],dtype=np.float64)*volume for f in fields])
bound=np.array([np.sum(f[...,6:8],dtype=np.float64)*volume for f in fields])
ax[1,1].plot(t,catalyst,'o-',color='#bf542b',label='catalyst present');ax[1,1].plot(t,bound,'o-',color='#227c85',label='bound food')
ax[1,1].set(title='Stored material',xlabel='time',ylabel='mass');ax[1,1].legend(frameon=False)
for x in ax.flat:x.spines[['top','right']].set_visible(False)
recipe_path=a.export.parent/'pilot.json'
recipe=json.loads(recipe_path.read_text())['recipe'] if recipe_path.exists() else {}
native=recipe.get('nvox','?')
fig.suptitle(f"Paid habitat pilot · native {native}³ · time {t[-1]:.3g}\n{last['coarse']}³ block means; {frames[0]['population']} founders",fontsize=14)
fig.savefig(a.output/'pilot-fields.png',dpi=170);fig.savefig(a.output/'pilot-fields.svg')
(a.output/'plotted-totals.json').write_text(json.dumps({'time':t.tolist(),'catalyst_present':catalyst.tolist(),'bound_food':bound.tolist(),'display':f"actual {last['coarse']}³ block means from native {native}³ GPU state; projections explicitly labeled"},indent=2)+'\n')
