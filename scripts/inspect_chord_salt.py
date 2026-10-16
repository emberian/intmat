"""Inspect actual 3D optical exports; compare matched-grid coupled/uncoupled controls."""
import argparse,json
from pathlib import Path
import numpy as np
from PIL import Image
p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('--control',type=Path);a=p.parse_args()
m=json.loads((a.run/'run.json').read_text());shape=tuple(m['shape'])
def field(root,name):return np.fromfile(root/'volume'/(name+'.f32'),dtype='<f4').reshape(shape,order='F')
h=field(a.run,'density');d=field(a.run,'deposit');dose=field(a.run,'activity');inside=h>.95
stats={'shape':shape,'mineral_quantiles_inside':np.quantile(d[inside],[0,.1,.25,.5,.75,.9,1]).tolist(),'fraction_support_above_isovalue':float(np.mean(d[inside]>.55)),'inventory_error':m['local_inventory_error']}
if a.control:
 c=field(a.control,'deposit');q=field(a.control,'activity');stats.update(mineral_mean_absolute_difference=float(np.mean(abs(c[inside]-d[inside]))),dose_relative_l2=float(np.linalg.norm((dose-q).ravel())/np.linalg.norm(dose.ravel())),mineral_correlation=float(np.corrcoef(c[inside],d[inside])[0,1]))
panels=[]
for arr in (d,dose/max(float(np.quantile(dose,.99)),1e-9)):
 for axis in range(3):
  z=np.clip(np.take(arr,shape[axis]//2,axis=axis).T,0,1);rgb=np.stack((z**.7,.6*z,z*.2),-1);panels.append(Image.fromarray((rgb*255).astype('uint8')).resize((320,320)))
im=Image.new('RGB',(960,640))
for i,v in enumerate(panels):im.paste(v,((i%3)*320,(i//3)*320))
im.save(a.run/'diagnostic.png');(a.run/'inspection.json').write_text(json.dumps(stats,indent=2));print(json.dumps(stats,indent=2))
