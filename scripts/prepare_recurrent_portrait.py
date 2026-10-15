"""Map a native recurrent material scalar into the existing volumetric renderer.
No simulation or invented microgeometry. Chosen ellipsoid is a viewing crop.
"""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('snapshot',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
a.out.mkdir(parents=True,exist_ok=False)
m=json.loads((a.snapshot/'state.json').read_text());shape=tuple(m['grid']['shape'])
def load(name):
 d=m['fields'][name];return np.fromfile(a.snapshot/d['file'],dtype=np.float32).reshape(shape,order='F')
deposit=load('deposit');activity=load('activity')
x,y,z=np.meshgrid(*[np.linspace(-1,1,n,dtype=np.float32) for n in shape],indexing='ij')
r=(x/.9)**2+(y/.48)**2+(z/.9)**2
mask=np.clip((1-r)/.04,0,1)
for name,field in [('structure',deposit),('activity',activity),('density',.02+.08*deposit),('mask',mask)]:
 field.astype('<f4').ravel(order='F').tofile(a.out/(name+'.f32'))
(a.out/'volume.json').write_text(json.dumps({'shape':shape,'bounds':[[-1,-1,-1],[1,1,1]],'fields':['structure','activity','density','mask'],'dtype':'float32','order':'F'},indent=2))
(a.out/'source-mapping.json').write_text(json.dumps({'source':str(a.snapshot.resolve()),'source_manifest_sha256':hashlib.sha256((a.snapshot/'state.json').read_bytes()).hexdigest(),'time':m['time'],'simulation_shape':shape,'structure':'unmodified recorded deposit','activity':'unmodified recorded fast recurrent activity','density':'.02+.08*deposit','viewing_crop':'ellipsoid axes .9,.48,.9; smooth4% radial edge; display selection only','deposit_quantiles':np.quantile(deposit,[0,.1,.5,.8,.9,1]).tolist()},indent=2))
