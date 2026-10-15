"""Native-resolution optical exports from dense solver fields; no resampling."""
import json, os, pathlib, shutil, sys
import numpy as np
root=pathlib.Path(sys.argv[1]); run=json.loads((root/'run.json').read_text()); volume=root/'volume';volume.mkdir(exist_ok=False)
surface_only=max(run['shape'])>=768
shape=run['shape']; count=int(np.prod(shape)); deposit=np.memmap(root/'deposit.f32',dtype='<f4',mode='r',shape=(count,))
for name in ('deposit','temperature','dose','pressure','intensity'):
    source=root/(name+'.f32')
    if source.exists(): os.link(source,volume/source.name)
for name in (() if surface_only else ('density','activity','structure','mask')):
    target=np.memmap(volume/(name+'.f32'),dtype='<f4',mode='w+',shape=(count,))
    for start in range(0,count,1048576):
        stop=min(start+1048576,count); dep=deposit[start:stop]
        target[start:stop]=.3+.45*dep if name=='density' else .025*dep if name=='activity' else dep if name=='structure' else 1.
    target.flush(); del target
fields=['deposit'] if surface_only else [p.stem for p in volume.glob('*.f32')]
(volume/'volume.json').write_text(json.dumps({'shape':shape,'bounds':[-1,-1,-1,1,1,1],'dtype':'little-endian float32','order':'F','fields':fields},indent=2))
job={'id':root.name,'family':'writing' if run['mode']=='writing' else 'waves','title':f'Native {shape[0]}³ thermal cursive' if run['mode']=='writing' else f'Native {shape[0]}³ acoustic mineral','description':f'Fixed extent 2³; actual native {shape[0]}³ simulation, no upsampling. Optical carrier is artistic; inclusion is measured deposit=.18.','volume':'volume','renderer':'volume','args':['--palette','amber-gold','--polished-carrier','--carrier-shape','ellipsoid','--inclusion-surface','deposit','--inclusion-isovalue','0.18','--emission','0.0','--absorption','0.16','--scattering','0.04','--surface-roughness','0.09','--orbit','0.4','--elevation','0.35'],'size':2400,'samples':128}
if surface_only: job['args']+=['--surface-only']
(root/'render-job.json').write_text(json.dumps(job,indent=2))
(root/'optical-export.json').write_text(json.dumps({'source':'deposit.f32','shape':shape,'resampling':False,'mapping':{'density':'.3+.45*deposit','activity':'.025*deposit','structure':'deposit','mask':'one; display carrier clipped by renderer'},'native_field_bytes':count*4},indent=2))
shutil.copyfile(run['config'],root/'recipe.toml')
print(json.dumps({'output':str(root),'shape':shape,'fields':fields,'render_job':job}))
