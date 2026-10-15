"""Render actual Julia electrical-growth states inside an artistic amber cutaway.
The carrier and optical transfers are staging; conductor topology and thermal
mineral/temperature fields come exclusively from saved simulation states.
"""
import argparse, gzip, hashlib, json, math, sys, time
from pathlib import Path
import bpy, numpy as np, openvdb as vdb
from mathutils import Vector
p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path)
p.add_argument('--steps',default='all');p.add_argument('--size',type=int,default=800);p.add_argument('--samples',type=int,default=48);p.add_argument('--device',default='METAL',choices=('METAL','CPU'));p.add_argument('--image-name',default='final.png');p.add_argument('--orbit',type=float,default=.38);p.add_argument('--elevation',type=float,default=.17);p.add_argument('--distance',type=float,default=4.4);p.add_argument('--halo-iso',type=float,default=.08);p.add_argument('--no-carrier',action='store_true');p.add_argument('--no-halo',action='store_true');p.add_argument('--halo-mode',choices=('surface','volume'),default='volume');p.add_argument('--carrier',choices=('cutaway','crystal'),default='crystal');p.add_argument('--backdrop',choices=('ink','twilight','slate'),default='ink');p.add_argument('--mineral-optics',choices=('clear','amber'),default='clear');p.add_argument('--save-scene',action='store_true')
a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);a.output.mkdir(parents=True,exist_ok=True)
renderer_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();run_sha=hashlib.sha256((a.run/'run.json').read_bytes()).hexdigest();run=json.loads((a.run/'run.json').read_text());steps=None if a.steps=='all' else set(map(int,a.steps.split(',')));frames=[f for f in run['frames'] if steps is None or f['step'] in steps]
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
scene=bpy.context.scene;scene.render.engine='CYCLES';scene.cycles.samples=a.samples;scene.cycles.use_denoising=True;scene.cycles.use_adaptive_sampling=True;scene.cycles.adaptive_threshold=.02;scene.cycles.adaptive_min_samples=min(a.samples,16);scene.cycles.max_bounces=16;scene.cycles.transmission_bounces=12;scene.cycles.volume_bounces=2
if a.device=='METAL':
 prefs=bpy.context.preferences.addons['cycles'].preferences;prefs.compute_device_type='METAL';prefs.get_devices()
 for d in prefs.devices:d.use=d.type=='METAL'
 assert any(d.use for d in prefs.devices);scene.cycles.device='GPU'
scene.render.resolution_x=a.size;scene.render.resolution_y=a.size;scene.render.resolution_percentage=100;scene.render.image_settings.file_format='PNG';scene.world.use_nodes=True;scene.world.node_tree.nodes['Background'].inputs['Color'].default_value=(.02,.025,.04,1);scene.world.node_tree.nodes['Background'].inputs['Strength'].default_value=.18
scene.view_settings.view_transform='AgX';scene.view_settings.look='AgX - Medium High Contrast';scene.view_settings.exposure=-.8

def aim(obj,at=(0,0,0)):obj.rotation_euler=(Vector(at)-obj.location).to_track_quat('-Z','Y').to_euler()
bpy.ops.object.camera_add(location=(a.distance*math.cos(a.elevation)*math.sin(a.orbit),-a.distance*math.cos(a.elevation)*math.cos(a.orbit),a.distance*math.sin(a.elevation)));cam=bpy.context.object;cam.data.lens=58;aim(cam);scene.camera=cam
for name,pos,power,color,size,sy in [('Honey softbox',(-3,-3,3),1000,(1,.66,.25),1.6,3),('Blue rim',(2,1,2),1100,(.24,.62,1),.8,3),('White ribbon',(-.5,1,3.5),1300,(1,.91,.76),.6,2.8),('Low warm',(0,-1,-2),80,(1,.19,.035),2,2)]:
 bpy.ops.object.light_add(type='AREA',location=pos);o=bpy.context.object;o.name=name;o.data.energy=power;o.data.color=color;o.data.shape='RECTANGLE';o.data.size=size;o.data.size_y=sy;aim(o)
def material(name,color,rough=.2,metal=0,trans=0,ior=1.45):
 m=bpy.data.materials.new(name);m.use_nodes=True;s=m.node_tree.nodes.get('Principled BSDF')
 for key,value in [('Base Color',color),('Roughness',rough),('Metallic',metal),('Transmission Weight',trans),('IOR',ior),('Coat Weight',.3)]:s.inputs[key].default_value=value
 return m,s
amber,ab=material('Cut amber optical carrier',(1,.87,.62,1),.055,0,1,1.47)
absorb=amber.node_tree.nodes.new('ShaderNodeVolumeAbsorption');absorb.inputs['Color'].default_value=(.72,.24,.035,1);absorb.inputs['Density'].default_value=.25;amber.node_tree.links.new(absorb.outputs[0],amber.node_tree.nodes.get('Material Output').inputs['Volume'])
halo,hs=material('Persistent thermally written mineral',(1,.83,.48,1),.07,0,1,1.14);hs.inputs['Thin Film Thickness'].default_value=285;hs.inputs['Thin Film IOR'].default_value=1.29
if a.mineral_optics=='amber':
 hs.inputs['IOR'].default_value=1.22;hs.inputs['Roughness'].default_value=.10
 absorption=halo.node_tree.nodes.new('ShaderNodeVolumeAbsorption');absorption.inputs['Color'].default_value=(.72,.24,.035,1);absorption.inputs['Density'].default_value=1.5;halo.node_tree.links.new(absorption.outputs[0],halo.node_tree.nodes.get('Material Output').inputs['Volume'])
copper,cs=material('Old conductive scars',(.62,.21,.035,1),.23,.78);cs.inputs['Emission Color'].default_value=(1,.2,.015,1);cs.inputs['Emission Strength'].default_value=.13
active,acs=material('Recently selected conductor sites',(.16,.65,1,1),.16,.55);acs.inputs['Emission Color'].default_value=(.035,.45,1,1);acs.inputs['Emission Strength'].default_value=2.5
seedmat,sm=material('Designed electrode roots',(.65,.3,.07,1),.18,.9)
black,bs=material('Black ceramic plinth',(.004,.006,.008,1),.3,.35)
# A faceted C-shaped cutaway, built from an extruded polygon. Its front window
# exposes the evolving interior. It is an optical observation chamber, not a
# simulated fracture or a predicted geometry.
polygon=[(-.91,-.62),(-.77,-.77),(-.57,-.77),(-.57,.37),(.58,.37),(.58,-.60),(.78,-.60),(.93,-.45),(.93,.54),(.77,.72),(-.76,.72),(-.91,.57)]
def prism(name,poly,z0,z1,mat):
 count=len(poly);vertices=[(x,y,z) for z in (z0,z1) for x,y in poly];faces=[tuple(reversed(range(count))),tuple(range(count,2*count))]+[(i,(i+1)%count,(i+1)%count+count,i+count) for i in range(count)]
 mesh=bpy.data.meshes.new(name);mesh.from_pydata(vertices,[],faces);mesh.update();o=bpy.data.objects.new(name,mesh);bpy.context.collection.objects.link(o);o.data.materials.append(mat)
 bevel=o.modifiers.new('Polished cut edge','BEVEL');bevel.width=.018;bevel.segments=3;o.modifiers.new('Weighted normals','WEIGHTED_NORMAL');return o
if not a.no_carrier:
 if a.carrier=='crystal':polygon=[(-.91,-.47),(-.64,-.70),(.48,-.64),(.87,-.26),(.82,.42),(.42,.68),(-.64,.62),(-.91,.26)]
 prism('Faceted amber observation chamber',polygon,-.89,.91,amber)
# prism('Thin floating basalt base',[(-.94,-.77),(.94,-.77),(.99,-.67),(.99,.72),(.88,.79),(-.88,.79),(-.99,.68),(-.99,-.67)],-.97,-.91,black)
bpy.ops.mesh.primitive_plane_add(size=200,location=(0,0,-1.005));floor_object=bpy.context.object;floor_object.data.materials.append(black)
if a.backdrop=='twilight':
 backdrop,bs2=material('Twilight blue photographic backdrop',(.055,.105,.16,1),.75)
 bpy.ops.mesh.primitive_plane_add(size=200,location=(0,2.2,0),rotation=(math.pi/2,0,0));bpy.context.object.data.materials.append(backdrop)
 scene.world.node_tree.nodes['Background'].inputs['Color'].default_value=(.22,.3,.45,1);scene.world.node_tree.nodes['Background'].inputs['Strength'].default_value=.35

if a.backdrop=='slate':
 bpy.data.objects.remove(floor_object,do_unlink=True)
 backdrop,bss=material('Matte blue photographic cove',(.025,.04,.06,1),1)
 bss.inputs['Specular IOR Level'].default_value=0;bss.inputs['Coat Weight'].default_value=0
 rows=[(-200,-1.005),(2,-1.005)]+[(2+3*math.sin(t),1.995-3*math.cos(t)) for t in np.linspace(0,math.pi/2,33)[1:]]+[(5,200)]
 mesh=bpy.data.meshes.new('Seamless slate cove');mesh.from_pydata([(x,y,z) for y,z in rows for x in (-200,200)],[],[(2*q,2*q+1,2*q+3,2*q+2) for q in range(len(rows)-1)]);mesh.update();obj=bpy.data.objects.new('Seamless slate cove',mesh);bpy.context.collection.objects.link(obj);mesh.materials.append(backdrop)
 for polygon in mesh.polygons:polygon.use_smooth=True
 scene.world.node_tree.nodes['Background'].inputs['Color'].default_value=(.1,.14,.2,1);scene.world.node_tree.nodes['Background'].inputs['Strength'].default_value=.35


def field_read(directory,name,shape):
 path=directory/(name+'.f32')
 if path.exists():return np.fromfile(path,dtype='<f4').reshape(shape,order='F')
 with gzip.open(str(path)+'.gz','rb') as stream:return np.frombuffer(stream.read(),dtype='<f4').reshape(shape,order='F')

def edge_paths(edges):
 adjacent={}
 for i,j in edges:
  adjacent.setdefault(int(i),[]).append(int(j));adjacent.setdefault(int(j),[]).append(int(i))
 unused={tuple(sorted((int(i),int(j)))) for i,j in edges};paths=[]
 def follow(start,next):
  path=[start];previous=start;current=next
  while True:
   unused.discard(tuple(sorted((previous,current))));path.append(current)
   if len(adjacent[current])!=2:break
   candidates=[j for j in adjacent[current] if tuple(sorted((current,j))) in unused]
   if not candidates:break
   previous,current=current,candidates[0]
  return path
 for start,neighbors in adjacent.items():
  if len(neighbors)!=2:
   for next in neighbors:
    if tuple(sorted((start,next))) in unused:paths.append(follow(start,next))
 while unused:
  start,following=iter(unused).__next__();paths.append(follow(start,following))
 return paths

def curves(name,positions,edges,radius,mat):
 data=bpy.data.curves.new(name,'CURVE');data.dimensions='3D';data.resolution_u=1;data.bevel_depth=radius;data.bevel_resolution=2
 for path in edge_paths(edges):
  points=positions[path]
  # Two corner-cutting passes only soften the lattice centerline by less than
  # one cell; endpoints/junctions and the entire graph connectivity are retained.
  for _ in range(2):
   if len(points)>2:
    q=.75*points[:-1]+.25*points[1:];r=.25*points[:-1]+.75*points[1:];points=np.vstack([points[0],np.column_stack([q,r]).reshape(-1,3),points[-1]])
  s=data.splines.new('POLY');s.points.add(len(points)-1)
  s.points.foreach_set('co',np.column_stack([points,np.ones(len(points))]).ravel())
 obj=bpy.data.objects.new(name,data);bpy.context.collection.objects.link(obj);data.materials.append(mat);obj['simulation_geometry']=True
 return obj

def mineral_volume(field,directory):
 grid=vdb.FloatGrid();grid.name='density';grid.copyFromArray(np.ascontiguousarray(field));grid.transform=vdb.createLinearTransform(voxelSize=2/field.shape[0]);grid.transform.postTranslate(tuple(-1+1/field.shape[0] for _ in range(3)))
 heat=vdb.FloatGrid();heat.name='temperature';heat.copyFromArray(np.ascontiguousarray(field_read(directory,'temperature',field.shape)));heat.transform=grid.transform
 path=a.output/f'mineral-{directory.name}.vdb';vdb.write(str(path),grids=[grid,heat]);data=bpy.data.volumes.new('Saved mineral density');data.filepath=str(path.resolve());obj=bpy.data.objects.new('Mineral memory volume',data);bpy.context.collection.objects.link(obj);obj['simulation_geometry']=True
 m=bpy.data.materials.get('Volumetric mineral memory')
 if m is None:
  m=bpy.data.materials.new('Volumetric mineral memory');m.use_nodes=True;nodes=m.node_tree.nodes;nodes.clear();output=nodes.new('ShaderNodeOutputMaterial');volume=nodes.new('ShaderNodeVolumePrincipled');volume.inputs['Color'].default_value=(1,.42,.07,1);volume.inputs['Density'].default_value=3.0;volume.inputs['Anisotropy'].default_value=.15;volume.inputs['Emission Color'].default_value=(.04,.5,1,1)
  attr=nodes.new('ShaderNodeAttribute');attr.attribute_name='density';multiply=nodes.new('ShaderNodeMath');multiply.operation='MULTIPLY';multiply.inputs[1].default_value=8;m.node_tree.links.new(attr.outputs['Fac'],multiply.inputs[0]);m.node_tree.links.new(multiply.outputs[0],volume.inputs['Density']);m.node_tree.links.new(volume.outputs['Volume'],output.inputs['Volume']);thermal=nodes.new('ShaderNodeAttribute');thermal.attribute_name='temperature';gain=nodes.new('ShaderNodeMath');gain.operation='MULTIPLY';gain.inputs[1].default_value=1.8;m.node_tree.links.new(thermal.outputs['Fac'],gain.inputs[0]);m.node_tree.links.new(gain.outputs[0],volume.inputs['Emission Strength'])
 data.materials.append(m);return {'grid':'density','mapping':'Saved mineral concentration times 8, amber scattering; saved temperature times 1.8 drives chosen blue emission' }

def isosurface(field,iso):
 grid=vdb.FloatGrid();grid.copyFromArray(np.ascontiguousarray(field));pts,tris,quads=grid.convertToPolygons(isovalue=iso,adaptivity=.04);pts=(pts+.5)*2/np.array(field.shape)-1
 mesh=bpy.data.meshes.new('Thermal mineral');mesh.from_pydata(pts.tolist(),[],tris.tolist()+quads.tolist());mesh.update();obj=bpy.data.objects.new('Thermal mineral',mesh);bpy.context.collection.objects.link(obj);obj.data.materials.append(halo);obj['simulation_geometry']=True
 for f in mesh.polygons:f.use_smooth=True
 return {'vertices':len(pts),'faces':len(tris)+len(quads),'isovalue':iso}
records=[]
for frame in frames:
 start=time.time();directory=a.run/frame['directory'];geometry=json.loads((directory/'geometry.json').read_text());shape=tuple(geometry['shape']);positions=(np.asarray(geometry['positions'],dtype=float)-.5)*2/np.array(shape)-1;edges=np.asarray(geometry['edges'],dtype=int)-1;birth=np.array(geometry['birth']);event=frame['step']
 for obj in list(bpy.data.objects):
  if obj.get('simulation_geometry'):
   data=obj.data;kind=obj.type;bpy.data.objects.remove(obj,do_unlink=True)
   (bpy.data.curves if kind=='CURVE' else bpy.data.volumes if kind=='VOLUME' else bpy.data.meshes).remove(data)
 if len(edges):
  age=event-np.maximum(birth[edges[:,0]],birth[edges[:,1]]);seed=(birth[edges[:,0]]==0)&(birth[edges[:,1]]==0);recent=(age<max(12,int(run['params']['capture_every']*1.5)))&~seed
  radius=max(.0015,.245*2/shape[0]);curves('Old electrical channels',positions,edges[~recent&~seed],radius,copper);curves('Recent electrical channels',positions,edges[recent],1.1*radius,active);curves('Designed electrode roots',positions,edges[seed],.009,seedmat)
 solid=field_read(directory,'solid',shape);surface_meta=(mineral_volume(solid,directory) if a.halo_mode=='volume' else isosurface(solid,a.halo_iso)) if not a.no_halo and solid.max()>a.halo_iso else None
 output=a.output/(a.image_name if len(frames)==1 else f"frame-{event:05d}.png");scene.render.filepath=str(output.resolve())
 input_files=[directory/'geometry.json']+[directory/(name+'.f32') if (directory/(name+'.f32')).exists() else directory/(name+'.f32.gz') for name in ('solid','temperature')]
 source_fields={str(path.relative_to(a.run)):hashlib.sha256(path.read_bytes()).hexdigest() for path in input_files}
 optical={'renderer_sha256':renderer_sha,'source':{'run_sha256':run_sha,'field_sha256':source_fields,'shape':shape,'step':event,'time':frame['time']},'settings':{**vars(a),'run':str(a.run),'output':str(a.output)},'camera':{'location':list(cam.location),'lens_mm':58,'fixed':True},'geometry':{'occupied_cells':len(positions),'edges':len(edges),'thermal_mineral':surface_meta},'interpretation':{'carrier':'No carrier.' if a.no_carrier else 'Designed faceted amber '+a.carrier+' is optical staging, not simulated fracture.','conductor':'Face-connected graph between occupied Julia cells; two subvoxel corner-cutting passes soften degree-two centerlines while preserving junctions and topology. Radius is a chosen visibility transfer. Cyan emission marks recently created edges, not physically solved luminescence.','mineral':('Isosurface of saved conservative thermal/electrical deposit, chosen translucent amber thin-film shader.' if a.halo_mode=='surface' else 'Saved mineral concentration times 8 drives volume scattering; saved temperature times 1.8 drives chosen blue emission.')+' No procedural animated deformation.'}}
 output.with_suffix('.png.optical.json').write_text(json.dumps(optical,indent=2))
 if a.save_scene:bpy.ops.wm.save_as_mainfile(filepath=str(output.with_suffix('.blend').resolve()))
 bpy.ops.render.render(write_still=True)
 record={'step':event,'time':frame['time'],'image':output.name,'sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'render_seconds':time.time()-start};records.append(record);(a.output/'frames.json').write_text(json.dumps({'source_run':str(a.run),'shape':shape,'fixed_camera':True,'frames':records},indent=2));print('FRAME_DONE',json.dumps(record),flush=True)
