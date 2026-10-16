"""Photograph actual coarse living patch state. Optical transfer only; no new dynamics."""
import argparse,json,hashlib,math,sys,time
from pathlib import Path
import bpy,numpy as np,openvdb as vdb
from mathutils import Vector
p=argparse.ArgumentParser();p.add_argument('report',type=Path);p.add_argument('output',type=Path);p.add_argument('--size',type=int,default=800);p.add_argument('--samples',type=int,default=32);p.add_argument('--time',type=float,default=64);p.add_argument('--device',default='METAL');p.add_argument('--density',type=float,default=40.0);p.add_argument('--trails',action='store_true');p.add_argument('--detection-floor',type=float,default=.12)
a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);a.output.parent.mkdir(parents=True,exist_ok=True);r=json.loads(a.report.read_text());f=min(r['frames'],key=lambda x:abs(x['time']-a.time));shape=r['metadata']['shape'];extent=np.array(r['metadata']['extent']);dx=extent/shape;shift=np.array([.75,-.5,0.]);cv=float(np.prod(dx));start=time.time()
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
s=bpy.context.scene;s.render.engine='CYCLES';s.cycles.samples=a.samples;s.cycles.use_denoising=True;s.cycles.volume_step_rate=.5;s.cycles.volume_bounces=6;s.cycles.max_bounces=10;s.render.resolution_x=s.render.resolution_y=a.size;s.render.resolution_percentage=100;s.render.image_settings.file_format='PNG';s.world.use_nodes=True;s.world.node_tree.nodes['Background'].inputs[0].default_value=(.035,.045,.065,1);s.world.node_tree.nodes['Background'].inputs[1].default_value=.55
if a.device=='METAL':
 prefs=bpy.context.preferences.addons['cycles'].preferences;prefs.compute_device_type='METAL';prefs.get_devices()
 for d in prefs.devices:d.use=d.type=='METAL'
 s.cycles.device='GPU'
s.view_settings.view_transform='AgX';s.view_settings.look='AgX - Medium High Contrast';s.view_settings.exposure=.2

def aim(o,at=(0,0,0)):o.rotation_euler=(Vector(at)-o.location).to_track_quat('-Z','Y').to_euler()
def mat(name,col,rough=.4,trans=0,metal=0):
 m=bpy.data.materials.new(name);m.use_nodes=True;b=m.node_tree.nodes.get('Principled BSDF');b.inputs['Base Color'].default_value=(*col,1);b.inputs['Roughness'].default_value=rough;b.inputs['Transmission Weight'].default_value=trans;b.inputs['Metallic'].default_value=metal;b.inputs['IOR'].default_value=1.35;return m
# Physical voxel centers use a consistent exact periodic translation.
fields={k:np.roll(np.array(v,dtype=np.float32).reshape(shape,order='F'),np.rint(shift/dx).astype(int),axis=(0,1,2))/cv for k,v in f['material'].items()}
# Single optical volume: total retained matrix sets density, A/B allocation sets color.
e=fields['catalyst_a']+fields['catalyst_b'];bound=fields['bound_a']+fields['bound_b'];density=e+bound*.25
ratio=(fields['catalyst_a']+.25*fields['bound_a'])/(density+1e-12)
density=np.maximum(0,density-a.detection_floor).astype(np.float32)
color=np.stack([.72-.60*ratio,.24+.28*ratio,.08+.55*ratio],axis=-1).astype(np.float32)
g=vdb.FloatGrid();g.name='density';g.copyFromArray(np.ascontiguousarray(density));g.transform=vdb.createLinearTransform(voxelSize=float(dx[0]));g.transform.postTranslate(tuple(-extent/2+dx/2))
c=vdb.Vec3SGrid();c.name='matrix_color';c.copyFromArray(np.ascontiguousarray(color));c.transform=g.transform
vp=a.output.with_suffix('.vdb');rg=vdb.FloatGrid();rg.name='reagent_ratio';rg.copyFromArray(np.ascontiguousarray(ratio.astype(np.float32)));rg.transform=g.transform;vdb.write(str(vp),grids=[g,c,rg]);bpy.ops.object.volume_add();vol=bpy.context.object;vol.name='Actual retained catalyst and bound-food field';vol.data.filepath=str(vp.resolve())
m=bpy.data.materials.new('Chosen optical absorption of retained chemistry');m.use_nodes=True;n=m.node_tree.nodes;n.clear();out=n.new('ShaderNodeOutputMaterial');pv=n.new('ShaderNodeVolumePrincipled');pv.inputs['Density'].default_value=a.density;pv.inputs['Density Attribute'].default_value='density';pv.inputs['Color Attribute'].default_value='matrix_color';pv.inputs['Color'].default_value=(1,1,1,1);pv.inputs['Anisotropy'].default_value=.25;da=n.new('ShaderNodeAttribute');da.attribute_name='density';mul=n.new('ShaderNodeMath');mul.operation='MULTIPLY';mul.inputs[1].default_value=a.density;m.node_tree.links.new(da.outputs['Fac'],mul.inputs[0]);m.node_tree.links.new(mul.outputs[0],pv.inputs['Density']);pv.inputs['Density Attribute'].default_value='';ca=n.new('ShaderNodeAttribute');ca.attribute_name='reagent_ratio';ramp=n.new('ShaderNodeValToRGB');ramp.color_ramp.elements[0].color=(.94,.76,.49,1);ramp.color_ramp.elements[1].color=(.46,.78,.94,1);m.node_tree.links.new(ca.outputs['Fac'],ramp.inputs[0]);m.node_tree.links.new(ramp.outputs[0],pv.inputs['Color']);em=n.new('ShaderNodeMath');em.operation='MULTIPLY';em.inputs[1].default_value=.85;m.node_tree.links.new(da.outputs['Fac'],em.inputs[0]);m.node_tree.links.new(em.outputs[0],pv.inputs['Emission Strength']);m.node_tree.links.new(ramp.outputs[0],pv.inputs['Emission Color']);pv.inputs['Color Attribute'].default_value='';m.node_tree.links.new(pv.outputs['Volume'],out.inputs['Volume']);vol.data.materials.append(m)
bodymats={}
for cell in f['cells']:
 energy=cell['energy']/cell['mass'];q=int(np.clip(energy,0,1)*15)
 if q not in bodymats:bodymats[q]=mat('Energy per mass optical bin '+str(q),(.82+.10*q/15,.78+.10*q/15,.68+.10*q/15),.2,1.0)
 bodymats[q].node_tree.nodes.get('Principled BSDF').inputs['IOR'].default_value=1.035
 bpy.ops.mesh.primitive_uv_sphere_add(segments=32,ring_count=20,radius=cell['radius'],location=(np.array(cell['position'])+shift)%extent-extent/2);o=bpy.context.object;o.name='Native body ID '+str(cell['id']);o.data.materials.append(bodymats[q])
 for face in o.data.polygons:face.use_smooth=True
bondmat=mat('Native mechanical bonds: chosen ivory optic',(.28,.17,.07),.45,.40,.05)
# Clip minimum-image bond segments to cube; replicate only the periodic segment copies that intersect it.
def clipped(p,q):
 d=q-p;lo=0.;hi=1.
 for k in range(3):
  if abs(d[k])<1e-12:
   if p[k]<0 or p[k]>extent[k]:return
  else:
   x=(0-p[k])/d[k];y=(extent[k]-p[k])/d[k];lo=max(lo,min(x,y));hi=min(hi,max(x,y))
 if hi<=lo:return
 return p+lo*d-extent/2,p+hi*d-extent/2
segments=0
for b in f['bonds']:
 p0=(np.array(b['a_position'])+shift)%extent;delta=np.array(b['b_position'])-np.array(b['a_position']);delta-=extent*np.round(delta/extent)
 for ix in (-1,0,1):
  for iy in (-1,0,1):
   for iz in (-1,0,1):
    pp=p0+extent*np.array([ix,iy,iz]);seg=clipped(pp,pp+delta)
    if seg is None:continue
    aa,bb=seg;diff=bb-aa;curve=bpy.data.curves.new('Native minimum-image junction','CURVE');curve.dimensions='3D';curve.bevel_depth=.008;curve.bevel_resolution=3;sp=curve.splines.new('POLY');sp.points.add(1);sp.points[0].co=(*aa,1);sp.points[1].co=(*bb,1);o=bpy.data.objects.new('Native bond '+str(b['a_id'])+'-'+str(b['b_id']),curve);bpy.context.collection.objects.link(o);curve.materials.append(bondmat);segments+=1
# Time-integrated observation of native cell motion: polygonal sampled paths, not deposited matter.
trailmat=mat('Observed cell paths, optical light writing',(.42,.24,.055),.35,.1,.3)
trailmat.node_tree.nodes.get('Principled BSDF').inputs['Emission Color'].default_value=(.25,.085,.012,1)
trailmat.node_tree.nodes.get('Principled BSDF').inputs['Emission Strength'].default_value=.25
trail_segments=0
history=[ff for ff in r['frames'] if f['time']-32<=ff['time']<=f['time']] if a.trails else []
for prev,nxt in zip(history,history[1:]):
 old={cc['id']:cc for cc in prev['cells']}
 for cc in nxt['cells']:
  if cc['id'] not in old:continue
  p0=(np.array(old[cc['id']]['position'])+shift)%extent;dd=np.array(cc['position'])-np.array(old[cc['id']]['position']);dd-=extent*np.round(dd/extent)
  for ix in (-1,0,1):
   for iy in (-1,0,1):
    for iz in (-1,0,1):
     pp=p0+extent*np.array([ix,iy,iz]);seg=clipped(pp,pp+dd)
     if seg is None:continue
     aa,bb=seg;curve=bpy.data.curves.new('Observed sampled trajectory','CURVE');curve.dimensions='3D';curve.bevel_depth=.0035;curve.bevel_resolution=2;sp=curve.splines.new('POLY');sp.points.add(1);sp.points[0].co=(*aa,1);sp.points[1].co=(*bb,1);o=bpy.data.objects.new('Past cell position ID '+str(cc['id']),curve);bpy.context.collection.objects.link(o);curve.materials.append(trailmat);trail_segments+=1
floor=mat('Midnight photographic floor',(.018,.025,.035),.8)
bpy.ops.mesh.primitive_plane_add(size=200,location=(0,0,-1.16));bpy.context.object.data.materials.append(floor)
for pos,power,col,size in [((-3,-4,5),500,(1,.73,.43),3),((3,1,4),750,(.5,.75,1),2),((0,3,1),500,(1,.35,.17),2)]:
 bpy.ops.object.light_add(type='AREA',location=pos);o=bpy.context.object;o.data.energy=power;o.data.color=col;o.data.shape='DISK';o.data.size=size;aim(o)
bpy.ops.object.camera_add(location=(2.3,-4.4,2.4));cam=bpy.context.object;cam.data.type='ORTHO';cam.data.ortho_scale=3.85;aim(cam);s.camera=cam
s.render.filepath=str(a.output.resolve());bpy.ops.render.render(write_still=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
meta=dict(title='After the last feeding',shape=shape,time=f['time'],cells=len(f['cells']),native_bonds=len(f['bonds']),rendered_periodic_segments=segments,observed_trajectory_segments=trail_segments,source_report_sha256=sha(a.report),source_frame_sha256=hashlib.sha256(json.dumps(f,sort_keys=True,separators=(',',':')).encode()).hexdigest(),renderer_sha256=sha(__file__),periodic_translation=shift.tolist(),voxel_spacing=dx.tolist(),voxel_volume=cv,optical_transfer={'retained_density':'catalyst_a+catalyst_b+.25*(bound_a+bound_b)','detection_floor':a.detection_floor,'scattering_scale':a.density,'fluorescence_scale':.85,'body_ior':1.035,'bond_radius':.008,'trajectory_overlay':a.trails},image_sha256=sha(a.output),settings={**vars(a),'report':str(a.report),'output':str(a.output)},render_seconds=time.time()-start,interpretation='Actual native-grid amounts divided by physical voxel volume; chosen fluorescence proportional to retained-field concentration, plus trilinear optical scattering of catalyst plus .25 bound food, after the recorded concentration detection floor, colored by A/B allocation. Native cell positions and physical body radii; body energy/mass selects optical color. Native ID-validated mechanical bonds shown at chosen optical radius .008; periodic minimum-image edges clipped to specimen bounds. Volume interpolation, opacity, color, bond thickness, light and floor are optical choices; Optional cell trails, when included, are piecewise-linear observations of actual sampled positions during the preceding32time units, not persistent deposited matter; exact periodic translation (.75,-.5,0) applies to all fields/bodies/bonds/trails. No additional simulated geometry, selected wavelength, adhesion, or higher simulation resolution is implied. No source fields are modified.')
a.output.with_suffix('.optical.json').write_text(json.dumps(meta,indent=2));print('DONE',json.dumps(meta),flush=True)
