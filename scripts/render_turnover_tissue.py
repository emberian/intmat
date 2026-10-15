"""Photograph native 3-D turnover fields; all optics/cutaways are display choices."""
import argparse, hashlib, json, math, sys, time
from pathlib import Path
import bpy
import numpy as np
import openvdb as vdb
from mathutils import Vector

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('run',type=Path);p.add_argument('output',type=Path)
p.add_argument('--steps',default='1500');p.add_argument('--size',type=int,default=720)
p.add_argument('--samples',type=int,default=48);p.add_argument('--device',choices=('CPU','METAL'),default='CPU')
p.add_argument('--bulk-iso',type=float,default=.25);p.add_argument('--matrix-iso',type=float,default=.16)
p.add_argument('--orbit',type=float,default=.62);p.add_argument('--elevation',type=float,default=.42)
p.add_argument('--distance',type=float,default=5.1);p.add_argument('--cut-x',type=float,default=-.12)
p.add_argument('--cut-y',type=float,default=.18);p.add_argument('--no-cutaway',action='store_true')
p.add_argument('--target',type=float,nargs=3,default=(0.,0.,0.))
p.add_argument('--image-name',default='final.png');p.add_argument('--save-scene',action='store_true')
a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);a.output.mkdir(parents=True,exist_ok=True)
run_bytes=(a.run/'run.json').read_bytes();run=json.loads(run_bytes);run_sha=hashlib.sha256(run_bytes).hexdigest()
renderer_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
steps=run['frame_steps'] if a.steps=='all' else list(map(int,a.steps.split(',')))
assert set(steps)<=set(run['frame_steps'])
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
scene=bpy.context.scene;scene.render.engine='CYCLES';scene.cycles.samples=a.samples
scene.cycles.use_denoising=True;scene.cycles.use_adaptive_sampling=True
scene.cycles.adaptive_threshold=.025;scene.cycles.adaptive_min_samples=min(16,a.samples)
scene.cycles.max_bounces=16;scene.cycles.transmission_bounces=12;scene.cycles.transparent_max_bounces=48
if a.device=='METAL':
 prefs=bpy.context.preferences.addons['cycles'].preferences;prefs.compute_device_type='METAL';prefs.get_devices()
 for d in prefs.devices:d.use=d.type=='METAL'
 assert any(d.use for d in prefs.devices);scene.cycles.device='GPU'
scene.render.resolution_x=scene.render.resolution_y=a.size;scene.render.resolution_percentage=100
scene.render.image_settings.file_format='PNG';scene.world.use_nodes=True
scene.world.node_tree.nodes['Background'].inputs['Color'].default_value=(.07,.08,.12,1)
scene.world.node_tree.nodes['Background'].inputs['Strength'].default_value=.3
scene.view_settings.view_transform='AgX';scene.view_settings.look='AgX - Medium High Contrast';scene.view_settings.exposure=-.45

def aim(o):o.rotation_euler=(-o.location).to_track_quat('-Z','Y').to_euler()
offset=(a.distance*math.cos(a.elevation)*math.sin(a.orbit),-a.distance*math.cos(a.elevation)*math.cos(a.orbit),a.distance*math.sin(a.elevation))
bpy.ops.object.camera_add(location=tuple(x+y for x,y in zip(a.target,offset)))
cam=bpy.context.object;cam.data.lens=55;cam.rotation_euler=(Vector(a.target)-cam.location).to_track_quat('-Z','Y').to_euler();scene.camera=cam
for name,pos,power,color,size in [('Ivory window',(-3,-4,4),1100,(1,.88,.70),3),('Lavender edge',(3,2,2.5),1000,(.52,.55,1),2),('Soft interior',(2,-3,.5),500,(.75,.9,1),2),('Vermilion rim',(-3,1,0),700,(1,.15,.07),2)]:
 bpy.ops.object.light_add(type='AREA',location=pos);o=bpy.context.object;o.name=name;o.data.energy=power;o.data.color=color;o.data.shape='DISK';o.data.size=size;aim(o)

def material(name,color,rough,metal,trans,ior):
 m=bpy.data.materials.new(name);m.use_nodes=True;s=m.node_tree.nodes.get('Principled BSDF')
 for k,v in [('Base Color',color),('Roughness',rough),('Metallic',metal),('Transmission Weight',trans),('IOR',ior)]:s.inputs[k].default_value=v
 return m,s

def clip_material(m,s):
 nodes,links=m.node_tree.nodes,m.node_tree.links;geom=nodes.new('ShaderNodeNewGeometry');xyz=nodes.new('ShaderNodeSeparateXYZ');links.new(geom.outputs['Position'],xyz.inputs[0]);conditions=[]
 for axis in 'XYZ':
  absolute=nodes.new('ShaderNodeMath');absolute.operation='ABSOLUTE';links.new(xyz.outputs[axis],absolute.inputs[0]);outside=nodes.new('ShaderNodeMath');outside.operation='GREATER_THAN';outside.inputs[1].default_value=1;links.new(absolute.outputs[0],outside.inputs[0]);conditions.append(outside.outputs[0])
 if not a.no_cutaway:
  x=nodes.new('ShaderNodeMath');x.operation='GREATER_THAN';x.inputs[1].default_value=a.cut_x;links.new(xyz.outputs['X'],x.inputs[0])
  y=nodes.new('ShaderNodeMath');y.operation='LESS_THAN';y.inputs[1].default_value=a.cut_y;links.new(xyz.outputs['Y'],y.inputs[0]);both=nodes.new('ShaderNodeMath');both.operation='MULTIPLY';links.new(x.outputs[0],both.inputs[0]);links.new(y.outputs[0],both.inputs[1]);conditions.append(both.outputs[0])
 condition=conditions[0]
 for other in conditions[1:]:
  maximum=nodes.new('ShaderNodeMath');maximum.operation='MAXIMUM';links.new(condition,maximum.inputs[0]);links.new(other,maximum.inputs[1]);condition=maximum.outputs[0]
 transparent=nodes.new('ShaderNodeBsdfTransparent');mix=nodes.new('ShaderNodeMixShader');links.new(condition,mix.inputs[0]);links.new(s.outputs[0],mix.inputs[1]);links.new(transparent.outputs[0],mix.inputs[2]);links.new(mix.outputs[0],nodes.get('Material Output').inputs['Surface'])

bulk,bs=material('Living plus persistent matter, ivory porcelain',(.72,.65,.48,1),.26,.08,.10,1.36)
attr=bulk.node_tree.nodes.new('ShaderNodeAttribute');attr.attribute_name='positive_phase'
bs.inputs['Emission Color'].default_value=(1,.055,.015,1)
bulk.node_tree.links.new(attr.outputs['Fac'],bs.inputs['Emission Strength']);clip_material(bulk,bs)
matrix,ms=material('Persistent paid matrix, violet glaze',(.12,.009,.055,1),.12,.28,.55,1.33)
ms.inputs['Thin Film Thickness'].default_value=420;ms.inputs['Thin Film IOR'].default_value=1.32;clip_material(matrix,ms)
back,bks=material('Charcoal cove',(.012,.014,.024,1),1,0,0,1.4);bks.inputs['Specular IOR Level'].default_value=0
rows=[(-200,-1.18),(3,-1.18)]+[(3+4*math.sin(t),2.82-4*math.cos(t)) for t in np.linspace(0,math.pi/2,33)[1:]]+[(7,200)]
mesh=bpy.data.meshes.new('Cove');mesh.from_pydata([(x,y,z) for y,z in rows for x in (-200,200)],[],[(2*q,2*q+1,2*q+3,2*q+2) for q in range(len(rows)-1)]);mesh.update();o=bpy.data.objects.new('Cove',mesh);bpy.context.collection.objects.link(o);mesh.materials.append(back)
for poly in mesh.polygons:poly.use_smooth=True

def surface(field,iso,material,name,excitation):
 grid=vdb.FloatGrid();grid.copyFromArray(np.ascontiguousarray(np.pad(field,1,mode='wrap')))
 points,tris,quads=grid.convertToPolygons(isovalue=iso,adaptivity=.015)
 indices=np.rint(points-1).astype(int)%np.array(field.shape)
 phase=np.maximum(0,excitation[tuple(indices.T)])*.45
 coords=(points-.5)*2/np.array(field.shape)-1
 mesh=bpy.data.meshes.new(name);mesh.from_pydata(coords.tolist(),[],tris.tolist()+quads.tolist());mesh.update()
 obj=bpy.data.objects.new(name,mesh);bpy.context.collection.objects.link(obj);obj['actual_field']=True;mesh.materials.append(material)
 for poly in mesh.polygons:poly.use_smooth=True
 attr=mesh.attributes.new('positive_phase','FLOAT','POINT');attr.data.foreach_set('value',phase.astype(np.float32))
 return {'vertices':len(points),'faces':len(tris)+len(quads),'threshold':iso}

records=[]
for step in steps:
 start=time.time();directory=a.run/'frames'/f'{step:07d}';meta=json.loads((directory/'volume.json').read_text());shape=tuple(meta['shape'])
 assert shape[2]>1 and meta['array_order']=='F' and meta['centering']=='cell'
 fields={};hashes={'volume.json':hashlib.sha256((directory/'volume.json').read_bytes()).hexdigest()}
 for name in ('biomass','matrix','excitation'):
  path=directory/(name+'.f32');raw=path.read_bytes();assert len(raw)==4*math.prod(shape);hashes[path.name]=hashlib.sha256(raw).hexdigest();fields[name]=np.frombuffer(raw,dtype='<f4').reshape(shape,order='F');assert np.isfinite(fields[name]).all()
 for obj in list(bpy.data.objects):
  if obj.get('actual_field'):
   data=obj.data;bpy.data.objects.remove(obj,do_unlink=True);bpy.data.meshes.remove(data)
 geometry={'bulk':surface(fields['biomass']+fields['matrix'],a.bulk_iso,bulk,'B plus M contour',fields['excitation']),'matrix':surface(fields['matrix'],a.matrix_iso,matrix,'M contour',fields['excitation'])}
 output=a.output/(a.image_name if len(steps)==1 else f'frame-{step:07d}.png');scene.render.filepath=str(output.resolve())
 optical={'renderer_sha256':renderer_sha,'source':{'run_sha256':run_sha,'step':step,'time':meta['time'],'shape':shape,'field_sha256':hashes},'settings':{**vars(a),'run':str(a.run),'output':str(a.output)},'camera':{'location':list(cam.location),'lens_mm':55,'target':list(a.target),'fixed':True},'geometry':geometry,'interpretation':'Two isosurfaces of native Julia concentration fields, B+M and M, with a fixed observation cutaway. Periodic padding is clipped at the physical box; no designed anatomy or extruded 2D structure. Ivory porcelain, violet glazing and weak positive-x emission are chosen optical mappings; no solved optics or luminescence. Thresholds, camera and phase gain remain fixed across time.'}
 output.with_suffix('.png.optical.json').write_text(json.dumps(optical,indent=2))
 if a.save_scene:bpy.ops.wm.save_as_mainfile(filepath=str(output.with_suffix('.blend').resolve()))
 bpy.ops.render.render(write_still=True)
 record={'step':step,'time':meta['time'],'image':output.name,'sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'render_seconds':time.time()-start};records.append(record);(a.output/'frames.json').write_text(json.dumps({'source_run_sha256':run_sha,'fixed_camera':True,'frames':records},indent=2));print('FRAME_DONE',json.dumps(record),flush=True)
