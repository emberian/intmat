#!/usr/bin/env python3
"""Plot retained native lattice observations; never advances a simulation."""
import argparse, csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--out',type=Path,default=Path('docs/figures'));a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
colors={'full':'#1874a5','no_mutation':'#d67920','no_transport_barrier':'#687778','no_recurrence':'#9852ab'}
labels={'full':'Full','no_mutation':'No mutation','no_transport_barrier':'No matrix barrier','no_recurrence':'No recurrence'}
fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained');rows=[]
fields=[('population','Living bodies'),('matrix','Matrix inventory'),('unique_genotypes','Living genotype codes'),('injury_region_population','Bodies in injury region'),('paid_work','Cumulative paid work'),('material_residual','Material ledger residual')]
for arm,c in colors.items():
 for seed in (11,29):
  data=[json.loads(x) for x in (a.root/f'seed-{seed}-{arm}'/'history.jsonl').read_text().splitlines()]
  for d in data:
   d['paid_work']=sum(d['work']); rows.append({'seed':seed,'arm':arm,**{k:d[k] for k in ['time']+[x[0] for x in fields]}})
  for ax,(key,title) in zip(axes.flat,fields):
   ax.plot([d['time'] for d in data],[d[key] for d in data],color=c,ls='-' if seed==11 else '--',lw=1.4,label=labels[arm] if seed==11 else None)
   ax.set(title=title,xlabel='Time');ax.axvline(80,color='black',lw=.6,alpha=.5);ax.grid(alpha=.15)
axes[0,0].legend(fontsize=8);fig.suptitle('Finite-budget lattice ecology: two seeds per condition\nSolid: seed 11; dashed: seed 29. Vertical line: conservative injury.',fontsize=12)
for ext in ('png','svg'):fig.savefig(a.out/f'lattice-ecology-trajectories.{ext}',dpi=160)
with (a.out/'lattice-ecology-trajectories.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
plt.close(fig)
run=a.root/'seed-11-full'/'frames';fig,axes=plt.subplots(4,3,figsize=(8,10),layout='constrained')
for r,(stem,title) in enumerate([('0000000-initial','Initial, t=0'),('0000800-pre-injury','Before injury, t=80'),('0000800-post-injury','After injury, t=80'),('0001600-sample','Final, t=160')]):
 meta=json.loads((run/f'{stem}.json').read_text());v=np.prod(meta['shape']);f=np.fromfile(run/f'{stem}.f32',dtype=np.float32).reshape((v,15),order='F');g=np.fromfile(run/f'{stem}.i32',dtype=np.int32).reshape((v,3),order='F')
 for c,(field,name,limit) in enumerate([(f[:,8],'Body mass',1.2),(f[:,3],'Matrix mass',1.2),(np.where(g[:,0]>0,(g[:,0]-1)%2,np.nan),'Inherited locus 1',1)]):
  ax=axes[r,c];im=ax.imshow(field.reshape((32,32),order='F').T,origin='lower',cmap='viridis' if c<2 else 'coolwarm',vmin=0,vmax=limit,interpolation='nearest');ax.set_xticks([]);ax.set_yticks([])
  if r==0:ax.set_title(name)
  if c==0:ax.set_ylabel(title)
  if r==3:fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.05,pad=.02)
fig.suptitle('Actual native state, seed 11 / full condition\nOne bit shows inherited code geography, not functional diversity.',fontsize=11)
for ext in ('png','svg'):fig.savefig(a.out/f'lattice-ecology-snapshots.{ext}',dpi=170)
