"""Inspect actual 3D material/current snapshots; projections are diagnostics only."""
import argparse,json
from pathlib import Path
import numpy as np
from PIL import Image,ImageDraw
p=argparse.ArgumentParser();p.add_argument('run',type=Path);a=p.parse_args();r=json.loads((a.run/'run.json').read_text());n=r['recipe']['n']
def field(step,name):return np.fromfile(a.run/'frames'/f'{step:05d}'/(name+'.f32'),dtype='<f4').reshape((n,n,n),order='F')
steps=[r['frames'][i]['step'] for i in [0,len(r['frames'])//2,3*len(r['frames'])//4,len(r['frames'])-1]]
fig=Image.new('RGB',(1280,1020),(16,21,30));draw=ImageDraw.Draw(fig)
colors={'deposit':[(12,20,28),(64,116,131),(230,232,220)],'activity':[(8,10,25),(117,31,111),(255,187,83)],'damage':[(6,20,31),(19,126,133),(213,238,99)]}
for col,step in enumerate(steps):
 for row,(name,vmax) in enumerate([('deposit',.8),('activity',.12),('damage',.7)]):
  arr=np.max(field(step,name),axis=1).T[::-1];u=np.clip(arr/vmax,0,1)*2;rgb=np.empty((*u.shape,3),dtype=np.uint8);stops=np.array(colors[name])
  for channel in range(3):rgb[:,:,channel]=np.interp(u,[0,1,2],stops[:,channel]).astype(np.uint8)
  im=Image.fromarray(rgb).resize((320,320));fig.paste(im,(col*320,row*340+20));draw.text((col*320+8,row*340+4),f'{name} · t={step*r["recipe"]["dt"]:.2f}',fill='white')
fig.save(a.run/'diagnostic.png')
late=field(steps[-1],'activity');prior=field(steps[-2],'activity');cut=np.quantile(prior,.98);hot=prior>=cut
report={'top2pct_previous_hotspot_final_retention':float(np.mean(late[hot]>=np.quantile(late,.98))),'previous_hotspot_mean_power_ratio':float(np.mean(late[hot])/np.mean(prior[hot])),'projection':'maximum along y; diagnostic only','steps':steps,'total_deposition':r['total_deposition'],'total_erosion':r['total_erosion'],'maximum_temperature':max(x['temperature_max'] for x in r['history']),'mass_residual':r['final_material']-r['initial_material'],'unconverged_solves':r['unconverged_solves']}
(a.run/'inspection.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
