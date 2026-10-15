#!/usr/bin/env python3
"""Apply declared optical staging to curated jobs, without touching simulation fields."""
import json,pathlib
root=pathlib.Path('outputs/collection-20261015')
# Only overrides selected optical settings; raw grids and run provenance remain untouched.
def configure(name,settings,unshelled=False):
 p=root/name/'render-job.json'
 if not p.exists():return
 j=json.loads(p.read_text());args=j['args']
 if unshelled and '--polished-carrier' in args:args.remove('--polished-carrier')
 for key,value in settings.items():
  if key in args:args[args.index(key)+1]=str(value)
  else:args.extend([key,str(value)])
 j['args']=args
 text=json.dumps(j,indent=2)+'\n'
 if p.read_text()!=text:p.write_text(text)
# Different camera angles, sections and materials make geometry legible rather than disguising it in a shared shell.
plans={
'phase-driven-coral-aged':({'--cutaway-axis':'y','--cutaway-position':'.35','--studio-background':'slate','--emission':0,'--density-scale':0,'--surface-smooth':1},False),
'phase-folded-strata-aged':({'--cutaway-axis':'x','--cutaway-position':'.25','--studio-background':'warm','--inclusion-tint':'amber','--density-scale':'.03','--emission':0,'--orbit':1.0,'--surface-smooth':1},False),
'phase-gyroid-cathedral-aged':({'--studio-background':'slate','--inclusion-tint':'jade','--density-scale':0,'--emission':0,'--surface-smooth':1,'--orbit':.25},True),
'phase-nested-interfaces-aged':({'--cutaway-axis':'y','--cutaway-position':0,'--studio-background':'warm','--inclusion-tint':'pearl','--density-scale':0,'--emission':0,'--surface-smooth':1,'--orbit':.1},False),
'phase-pearl-orchard-aged':({'--studio-background':'slate','--ior':1.18,'--density-scale':.03,'--emission':0,'--surface-smooth':1,'--inclusion-tint':'rose'},False),
'phase-perforated-veil-open':({'--studio-background':'warm','--density-scale':0,'--emission':0,'--surface-smooth':1,'--inclusion-tint':'jade','--orbit':.8},True),
'write-counterspiral':({'--studio-background':'warm','--density-scale':0,'--elevation':.2},True),
'write-engulfed-script':({'--studio-background':'warm','--absorption':.05,'--scattering':.01,'--ior':1.25,'--inclusion-tint':'rose'},False),
'write-nested-script':({'--cutaway-axis':'y','--cutaway-position':.4,'--studio-background':'slate','--density-scale':0,'--inclusion-tint':'amber'},False),
'write-pearl-nuclei':({'--studio-background':'slate','--ior':1.25,'--absorption':.035,'--scattering':.01},False),
'write-pulse-rosary':({'--studio-background':'warm','--density-scale':0,'--inclusion-tint':'rose'},True),
'write-trefoil':({'--studio-background':'slate','--density-scale':0,'--inclusion-tint':'jade','--orbit':.2},True),
'branch-thunder-root':({'--density-scale':.12,'--emission':.3,'--absorption':.12,'--scattering':.025,'--camera-distance':4.9},False),
'branch-copper-thicket':({'--density-scale':0,'--emission':.15,'--studio-background':'slate','--graph-radius':.004,'--orbit':.3},True),
'branch-fault-script':({'--density-scale':.05,'--emission':.15,'--z-scale':.5,'--studio-background':'warm','--ior':1.2,'--elevation':.65},False),
'branch-field-comb':({'--density-scale':0,'--emission':.15,'--studio-background':'warm','--graph-radius':.004,'--orbit':1.1},True),
'branch-ice-antlers':({'--density-scale':0,'--emission':.15,'--studio-background':'slate','--graph-radius':.004,'--orbit':.7},True),
'branch-three-veins':({'--density-scale':.04,'--emission':.15,'--studio-background':'slate','--ior':1.2,'--camera-distance':4.9},False),
}
for name,(settings,unshelled) in plans.items():configure(name,settings,unshelled)
mechanics=['golden-sail','ice-cleavage','shear-veil','copper-lace','parted-choir','held-breath']
for i,name in enumerate(mechanics):
 configure('fracture-'+name,{'--camera-distance':3.2 if i in(0,2,3,5) else 3.8,'--elevation':.35 if i==4 else (.55 if i%2==0 else .32),'--studio-background':'warm' if i%2==0 else 'slate','--density-scale':.12,'--orbit':.1 if i==4 else .4+i*.17})
# Use the measured pressure crest scalar as a surface instead of a diffuse glowing orb.
waves=['crossed-choir','memory-rings','neural-lattice','neural-scroll','oblique-fronts','pearl-cavity']
for i,name in enumerate(waves):
 configure('waves-'+name,{'--inclusion-surface':'activity','--inclusion-isovalue':[.38,.45,.3,.4,.4,.22][i],'--surface-smooth':1,'--density-scale':0,'--emission':0,'--studio-background':'slate' if i%2==0 else 'warm','--camera-distance':[3.0,4.0,4.0,3.8,4.0,4.2][i],'--inclusion-tint':['pearl','amber','jade','rose','pearl','amber'][i],'--orbit':.2+i*.18,'--cutaway-axis':'y' if i%2==0 else 'x','--cutaway-position':[.3,0,.4,.2,.1,1.0][i]},True)

configure('ecology-tidal-lamellae',{'--density-scale':.04,'--emission':.02,'--studio-background':'warm','--z-scale':.8,'--inclusion-surface':'phase','--inclusion-isovalue':0,'--inclusion-transmission':.88,'--inclusion-tint':'jade','--elevation':.2},True)

for selection in root.glob('collection-*-selected.json'):
 data=json.loads(selection.read_text());ids=data if isinstance(data,list) else data.get('ids',[])
 for ident in ids:
  p=root/ident/'render-job.json'
  if p.exists() and not ident.startswith('ecology-'):
   j=json.loads(p.read_text())
   if '--frame-content' not in j['args']:
    j['args'].append('--frame-content');p.write_text(json.dumps(j,indent=2)+'\n')
configure('ecology-pearl-niches',{'--density-scale':.02,'--emission':.03,'--studio-background':'slate','--inclusion-surface':'phase','--inclusion-isovalue':.2,'--inclusion-transmission':.65,'--inclusion-tint':'pearl','--surface-smooth':1,'--orbit':.1},True)
configure('ecology-root-canals',{'--density-scale':.12,'--emission':.03,'--studio-background':'warm','--z-scale':.45,'--ior':1.22,'--elevation':.65})
configure('ecology-finite-orchard',{'--density-scale':.07,'--emission':.02,'--studio-background':'slate','--carrier-shape':'ellipsoid','--ior':1.2,'--orbit':1.1})
configure('ecology-crossfeed-foam',{'--density-scale':.015,'--emission':.02,'--studio-background':'warm','--inclusion-surface':'phase','--inclusion-isovalue':0,'--inclusion-transmission':.6,'--inclusion-tint':'rose','--surface-smooth':1,'--orbit':.3},True)
