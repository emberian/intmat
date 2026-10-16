"""Photograph native finite-precipitation geometry; no designed interior surface."""
import argparse,hashlib,json,math,sys,time
from pathlib import Path
import bpy,numpy as np,openvdb as vdb
from mathutils import Vector
p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path);p.add_argument('--size',type=int,default=720);p.add_argument('--samples',type=int,default=32);p.add_argument('--device',choices=('CPU','METAL'),default='CPU');p.add_argument('--step',type=int);p.add_argument('--mineral-iso',type=float,default=.45);p.add_argument('--wet-iso',type=float,default=.10);p.add_argument('--no-wet',action='store_true');p.add_argument('--orbit',type=float,default=.45);p.add_argument('--elevation',type=float,default=.63)
a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);a.output.parent.mkdir(parents=True,exist_ok=True);renderer_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();r=json.loads((a.run/'run.json').read_text());frame=r['frames'][-1] if a.step is None else next(f for f in r['frames'] if f['step']==a.step);directory=a.run/frame['directory'];shape=tuple(r['shape'])
def read(name):return np.fromfile(directory/(name+'.f32'),dtype='<f4').reshape(shape,order='F')
def sample(field,points):
 q=np.clip((points+1)*np.array(shape)/2-.5,0,np.array(shape)-1);lo=np.floor(q).astype(int);hi=np.minimum(lo+1,np.array(shape)-1);weight=q-lo;result=np.zeros(len(points))
 for i in (0,1):
  for j in (0,1):
   for k in (0,1):
    indices=[(hi if bit else lo)[:,d] for d,bit in enumerate((i,j,k))];w=np.prod(np.column_stack([weight[:,d] if bit else 1-weight[:,d] for d,bit in enumerate((i,j,k))]),axis=1);result+=field[tuple(indices)]*w
 return result
mineral=read('mineral');species_a=read('a');species_b=read('b');phi=read('phi')
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
s=bpy.context.scene;s.render.engine='CYCLES';s.cycles.samples=a.samples;s.cycles.use_denoising=True;s.cycles.max_bounces=12;s.cycles.transmission_bounces=8;s.cycles.use_adaptive_sampling=True;s.cycles.adaptive_threshold=.02
if a.device=='METAL':
 prefs=bpy.context.preferences.addons['cycles'].preferences;prefs.compute_device_type='METAL';prefs.get_devices()
 for d in prefs.devices:d.use=d.type=='METAL'
 s.cycles.device='GPU'
s.render.resolution_x=a.size;s.render.resolution_y=a.size;s.render.resolution_percentage=100;s.render.image_settings.file_format='PNG';s.world.use_nodes=True;s.world.node_tree.nodes['Background'].inputs['Color'].default_value=(.35,.39,.34,1);s.world.node_tree.nodes['Background'].inputs['Strength'].default_value=.28;s.view_settings.view_transform='AgX';s.view_settings.look='AgX - Medium High Contrast';s.view_settings.exposure=-.7

def aim(o,at=(0,0,0)):o.rotation_euler=(Vector(at)-o.location).to_track_quat('-Z','Y').to_euler()
bpy.ops.object.camera_add(location=(4*math.cos(a.elevation)*math.sin(a.orbit),-4*math.cos(a.elevation)*math.cos(a.orbit),4*math.sin(a.elevation)));cam=bpy.context.object;cam.data.lens=62;aim(cam);s.camera=cam
for pos,power,color,size,sy in [((-2,-3,4),500,(1,.91,.75),2.5,3.5),((2,1,2),650,(.61,.82,1),1.2,3),((0,2,4),800,(1,1,.94),.5,2.5)]:
 bpy.ops.object.light_add(type='AREA',location=pos);o=bpy.context.object;o.data.energy=power;o.data.color=color;o.data.shape='RECTANGLE';o.data.size=size;o.data.size_y=sy;aim(o)
def material(name,color,roughness,metallic=0,transmission=0,ior=1.45):
 m=bpy.data.materials.new(name);m.use_nodes=True;bs=m.node_tree.nodes.get('Principled BSDF')
 for key,value in [('Base Color',color),('Roughness',roughness),('Metallic',metallic),('Transmission Weight',transmission),('IOR',ior),('Coat Weight',.3)]:bs.inputs[key].default_value=value
 return m,bs
copper,cb=material('Finite mineral, residual-reagent color',(.08,.28,.19,1),.32,.35)
attribute=copper.node_tree.nodes.new('ShaderNodeAttribute');attribute.attribute_name='reagent_optics';copper.node_tree.links.new(attribute.outputs['Color'],cb.inputs['Base Color']);cb.inputs['Coat Weight'].default_value=.04;cb.inputs['Thin Film Thickness'].default_value=0;cb.inputs['Thin Film IOR'].default_value=1.38
wet,wb=material('Low concentration wet interface',(.7,.94,.83,1),.045,0,1,1.085)
back,bb=material('Quiet sage photographic cove',(.26,.29,.24,1),.95);bb.inputs['Coat Weight'].default_value=0;bb.inputs['Specular IOR Level'].default_value=.1
rough=copper.node_tree.nodes.new('ShaderNodeAttribute');rough.attribute_name='phase_roughness';copper.node_tree.links.new(rough.outputs['Fac'],cb.inputs['Roughness'])
metal=copper.node_tree.nodes.new('ShaderNodeAttribute');metal.attribute_name='residual_metal';copper.node_tree.links.new(metal.outputs['Fac'],cb.inputs['Metallic'])
allpoints=[];geometry=[]
surfaces=[(a.mineral_iso,'Saved precipitated mineral',copper)]
if not a.no_wet:surfaces.insert(0,(a.wet_iso,'Saved low concentration wet interface',wet))
for iso,name,mat in surfaces:
 grid=vdb.FloatGrid();grid.copyFromArray(np.ascontiguousarray(mineral));pts,tris,quads=grid.convertToPolygons(isovalue=iso,adaptivity=.015);pts=(pts+.5)*2/np.array(shape)-1
 mesh=bpy.data.meshes.new(name);mesh.from_pydata(pts.tolist(),[],tris.tolist()+quads.tolist());mesh.update();o=bpy.data.objects.new(name,mesh);bpy.context.collection.objects.link(o);mesh.materials.append(mat)
 for face in mesh.polygons:face.use_smooth=True
 if mat==copper:
  av=sample(species_a,pts);bv=sample(species_b,pts);ratio=av/(av+bv+1e-8);colors=np.column_stack([.10+.445*ratio,.25-.055*ratio,.19-.12*ratio,np.ones(len(pts))]);attr=mesh.color_attributes.new(name='reagent_optics',type='FLOAT_COLOR',domain='POINT');attr.data.foreach_set('color',colors.ravel());phase_abs=np.minimum(1,np.abs(sample(phi,pts)));rough_attr=mesh.attributes.new('phase_roughness','FLOAT','POINT');rough_attr.data.foreach_set('value',.40+.22*phase_abs);metal_attr=mesh.attributes.new('residual_metal','FLOAT','POINT');metal_attr.data.foreach_set('value',.12+.50*ratio)
 allpoints.append(pts);geometry.append({'native_mineral_isovalue':iso,'vertices':len(pts),'faces':len(tris)+len(quads)})
points=np.vstack(allpoints);minimum_z=float(np.min(points[:,2]));bpy.ops.mesh.primitive_plane_add(size=200,location=(0,0,minimum_z-.018));bpy.context.object.data.materials.append(back)
s.render.filepath=str(a.output.resolve());start=time.time();bpy.ops.render.render(write_still=True)
optical={'title':r['title'],'shape':shape,'step':frame['step'],'time':frame['time'],'renderer_sha256':renderer_sha,'source_run_sha256':hashlib.sha256((a.run/'run.json').read_bytes()).hexdigest(),'source_fields_sha256':{name:hashlib.sha256((directory/(name+'.f32')).read_bytes()).hexdigest() for name in ('mineral','a','b','phi')},'image_sha256':hashlib.sha256(a.output.read_bytes()).hexdigest(),'geometry':geometry,'settings':{**vars(a),'run':str(a.run),'output':str(a.output)},'render_seconds':time.time()-start,'interpretation':'All displayed surfaces are native mineral concentration isosurfaces; no render mask, surface deformation or invented internal mesh. Coordinates are uniformly divided by physical half-extent for photographic staging. Residual A/(A+B) trilinearly sampled at native mineral vertices selects copper-to-verdigris color; saved absolute phase sets roughness 0.40–0.62 and residual ratio sets metallic response 0.12–0.62. These are chosen optical mappings, not an oxidation solve. The selected surface shaders and photographic floor/lights are optical choices, not calibrated material optics. A lower-concentration clear film is included only when no_wet is false.'}
a.output.with_suffix('.optical.json').write_text(json.dumps(optical,indent=2));print('DONE',json.dumps({'render_seconds':optical['render_seconds'],'geometry':geometry}),flush=True)
