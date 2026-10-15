#!/usr/bin/env python3
"""Review actual native 3D lattice frames and export unsmoothed scalar volumes."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path);p.add_argument('--injury-step',type=int,default=1600);p.add_argument('--final-step',type=int,default=3200);p.add_argument('--export-fields',action='store_true');a=p.parse_args()
if a.output.exists():raise SystemExit('fresh review directory required')
a.output.mkdir(parents=True);exports=[];metrics=[]
states=[(0,'initial'),(a.injury_step,'pre-injury'),(a.injury_step,'post-injury'),(a.final_step,'sample')]
fig,axes=plt.subplots(4,3,figsize=(8.5,11),layout='constrained')
for row,(step,label) in enumerate(states):
 stem=f'{step:07d}-{label}';frame=a.run/'frames'/stem;m=json.loads(frame.with_suffix('.json').read_text());dims=tuple(m['shape']);v=np.prod(dims)
 f=np.fromfile(frame.with_suffix('.f32'),np.float32).reshape((v,15),order='F');ints=np.fromfile(frame.with_suffix('.i32'),np.int32).reshape((v,3),order='F');live=ints[:,0]>0
 b=f[:,8].reshape(dims,order='F');matrix=f[:,3].reshape(dims,order='F');g=ints[:,0].reshape(dims,order='F')
 neighbor_total=neighbor_same=0
 for axis in range(3):
  other=np.roll(g,1,axis=axis);both=(g>0)&(other>0);neighbor_total+=int(both.sum());neighbor_same+=int(((g==other)&both).sum())
 codes,counts=np.unique(ints[live,0],return_counts=True);mixed=float(np.sum((counts/max(1,live.sum()))**2))
 metrics.append(dict(step=step,label=label,time=m['time'],population=int(live.sum()),living_genotype_codes=len(codes),
  neighbor_live_edges=neighbor_total,neighbor_same_code_fraction=neighbor_same/max(1,neighbor_total),
  abundance_matched_independent_code_agreement=mixed,interpretation='Local inherited-code geography; not functional phenotype or cooperative individuality.'))
 for col,(field,title,vmax,cmap) in enumerate([(b,'Body mass',1.2,'viridis'),(matrix,'Matrix mass',1.2,'viridis'),(np.where(g>0,(g-1)%2,np.nan),'Inherited locus 1',1,'coolwarm')]):
  ax=axes[row,col];im=ax.imshow(field[:,:,dims[2]//2].T,origin='lower',interpolation='nearest',vmin=0,vmax=vmax,cmap=cmap)
  ax.set_xticks([]);ax.set_yticks([])
  if row==0:ax.set_title(title)
  if col==0:ax.set_ylabel(f't={m["time"]:.0f}, {label}')
  if row==3:fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.05,pad=.02)
 for index,name in [(8,'biomass'),(3,'matrix'),(2,'precursor'),(9,'stored_work')]:
  path=a.output/f'{stem}-{name}.f32';field=np.asarray(f[:,index],dtype='<f4')
  if a.export_fields:path.write_bytes(field.tobytes())
  exports.append(dict(path=path.name if a.export_fields else None,source_channel_index_zero_based=index,field=name,dimensions=dims,array_order='F',dtype='little-endian Float32',extent=dims,spacing=[1,1,1],
   source=str(frame.with_suffix('.f32').resolve()),source_sha256=hashlib.sha256(frame.with_suffix('.f32').read_bytes()).hexdigest(),
   sha256=hashlib.sha256(field.tobytes()).hexdigest(),minimum=float(field.min()),maximum=float(field.max()),step=step,time=m['time'],label=label))
fig.suptitle('Actual 32³ lattice volume: central z section, no interpolation\nSpherical injury is imposed; one inherited bit is a diagnostic code label.',fontsize=11)
for ext in ('png','svg'):fig.savefig(a.output/f'native-sections.{ext}',dpi=180)
captures=[]
for meta in sorted((a.run/'frames').glob('*.json')):
 m=json.loads(meta.read_text());files={}
 for ext in ('.json','.f32','.i32'):
  path=meta.with_suffix(ext);files[ext]=dict(path=str(path.resolve()),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
 captures.append(dict(metadata=m,files=files))
(a.output/'field-contract.json').write_text(json.dumps(dict(captures=captures,native_dimensions=[32,32,32],exports=exports,observations=metrics,
 interpretation='Scalar exports preserve exact native values and spacing. No optical transfer, contour smoothing, invented geometry or resolution enhancement applied.'),indent=2)+'\n')
print(json.dumps(metrics,indent=2))
