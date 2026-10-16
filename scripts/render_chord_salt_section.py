"""Photograph actual wave-written mineral inside its authored finite precursor support."""
import argparse,json,sys,hashlib,math
from pathlib import Path
import bpy,bmesh,numpy as np,openvdb as vdb
from mathutils import Vector
p=argparse.ArgumentParser();p.add_argument('frame',type=Path);p.add_argument('output',type=Path);p.add_argument('--no-floor',action='store_true');p.add_argument('--size',type=int,default=640);p.add_argument('--samples',type=int,default=16);p.add_argument('--iso',type=float,default=.55);p.add_argument('--device',choices=('CPU','METAL'),default='CPU');p.add_argument('--threads',type=int,default=2);a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
if a.threads<1:p.error('--threads must be positive')
meta=json.loads((a.frame/'volume.json').read_text());shape=tuple(meta['shape']);nodal='sample_convention' in meta;h=2/(shape[0]-1) if nodal else 2/shape[0]
fields={name:np.fromfile(a.frame/(name+'.f32'),dtype='<f4').reshape(shape,order='F') for name in ['deposit','activity','density']}
grid=vdb.FloatGrid();grid.copyFromArray(np.ascontiguousarray(fields['deposit']));points,tri,quad=grid.convertToPolygons(isovalue=a.iso,adaptivity=.05)
indices=np.clip(np.rint(points).astype(int),0,np.array(shape)-1);points=points*h-1 if nodal else (points+.5)*h-1
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
mesh=bpy.data.meshes.new('Computed mineral level set');faces=tri.tolist()+quad.tolist();mesh.from_pydata(points.tolist(),[],faces);mesh.update();obj=bpy.data.objects.new('Conductive mineral / actual3D state',mesh);bpy.context.collection.objects.link(obj)
bm=bmesh.new();bm.from_mesh(mesh);bmesh.ops.bisect_plane(bm,geom=list(bm.verts)+list(bm.edges)+list(bm.faces),dist=1e-6,plane_co=(0,-.04,0),plane_no=(0,1,0),clear_inner=True,clear_outer=False);bm.to_mesh(mesh);bm.free()
for poly in mesh.polygons:poly.use_smooth=True
def mat(name,color,rough,trans,ior):
 m=bpy.data.materials.new(name);m.use_nodes=True;s=m.node_tree.nodes.get('Principled BSDF');s.inputs['Base Color'].default_value=color;s.inputs['Roughness'].default_value=rough;s.inputs['Transmission Weight'].default_value=trans;s.inputs['IOR'].default_value=ior;s.inputs['Coat Weight'].default_value=.1;return m,s
m,s=mat('Opaline conductor / current and damage optics',(.065,.16,.22,1),.34,.08,1.48);s.inputs['Metallic'].default_value=.08;s.inputs['Thin Film IOR'].default_value=1.33;s.inputs['Emission Color'].default_value=(.08,.27,.5,1)
mesh.materials.append(m)
# The translucent precursor surface is the actual authored finite support, not a box.
host=vdb.FloatGrid();host.copyFromArray(np.ascontiguousarray(fields['density']));hp,ht,hq=host.convertToPolygons(isovalue=.45,adaptivity=.01);hp=(hp+.5)*h-1
me=bpy.data.meshes.new('Finite precursor / actual support');hfaces=ht.tolist()+hq.tolist();me.from_pydata(hp.tolist(),[],hfaces);me.update();bm=bmesh.new();bm.from_mesh(me);bmesh.ops.bisect_plane(bm,geom=list(bm.verts)+list(bm.edges)+list(bm.faces),dist=1e-6,plane_co=(0,.12,0),plane_no=(0,1,0),clear_inner=True,clear_outer=False);bm.to_mesh(me);bm.free();o=bpy.data.objects.new('Translucent salt precursor',me);bpy.context.collection.objects.link(o)
for face in me.polygons:face.use_smooth=True
medium,ms=mat('Honey optical precursor',(.89,.55,.16,1),.055,1,1.37);absorb=medium.node_tree.nodes.new('ShaderNodeVolumeAbsorption');absorb.inputs['Color'].default_value=(.75,.31,.075,1);absorb.inputs['Density'].default_value=.32;medium.node_tree.links.new(absorb.outputs[0],medium.node_tree.nodes.get('Material Output').inputs['Volume']);me.materials.append(medium)
sc=bpy.context.scene;sc.render.engine='CYCLES';sc.cycles.device='CPU';sc.render.threads_mode='FIXED';sc.render.threads=a.threads;
if a.device=='METAL':
 prefs=bpy.context.preferences.addons['cycles'].preferences;prefs.compute_device_type='METAL';prefs.get_devices()
 for device in prefs.devices:device.use=device.type=='METAL'
 if not any(device.use for device in prefs.devices):raise RuntimeError('No enabled Metal device')
 sc.cycles.device='GPU'
sc.cycles.samples=a.samples;sc.cycles.use_denoising=True;sc.cycles.max_bounces=16;sc.cycles.transmission_bounces=12
sc.render.resolution_x=a.size;sc.render.resolution_y=a.size;sc.render.resolution_percentage=100;sc.world.use_nodes=True;sc.world.node_tree.nodes['Background'].inputs[0].default_value=(.10,.16,.23,1);sc.world.node_tree.nodes['Background'].inputs[1].default_value=.3;sc.view_settings.view_transform='AgX';sc.view_settings.look='AgX - Medium High Contrast';sc.view_settings.exposure=.1
def aim(o,at=(0,0,0)):o.rotation_euler=(Vector(at)-o.location).to_track_quat('-Z','Y').to_euler()
bpy.ops.object.camera_add(location=(2.7,-4,2.1));camera=bpy.context.object;camera.data.type='ORTHO';camera.data.ortho_scale=2.25;aim(camera);sc.camera=camera
for name,pos,power,color,size,sy in [('Warm silk',(-3,-3,4),650,(1,.88,.73),3,3),('White edge',(3,1,3),700,(.76,.9,1),2,3),('Cool fill',(-2,1,.2),170,(.6,.8,1),2,2.8)]:
 bpy.ops.object.light_add(type='AREA',location=pos);o=bpy.context.object;o.name=name;o.data.energy=power;o.data.color=color;o.data.shape='RECTANGLE';o.data.size=size;o.data.size_y=sy;aim(o)
floor,fs=mat('Slate photographic cove',(.06,.075,.085,1),.95,0,1.45);bpy.ops.mesh.primitive_plane_add(size=200,location=(0,0,-.88));bpy.context.object.data.materials.append(floor);bpy.context.object.hide_render=a.no_floor
a.output.parent.mkdir(parents=True,exist_ok=True);sc.render.image_settings.file_format='PNG';sc.render.filepath=str(a.output.resolve());bpy.ops.render.render(write_still=True)
a.output.with_suffix('.optical.json').write_text(json.dumps({'source_frame':a.frame.name,'fields':'Actual wave-written mineral levelset, exactly plane-clipped at y=-.04 without caps. Assigned matte navy optics; no field smoothing or emission.','source_shape':shape,'node_centered':nodal,'isovalue':a.iso,'vertices':len(points),'polygons':len(mesh.polygons),'cut':'Exact bmesh plane clipping, open without caps','floor':not a.no_floor,'carrier':'Actual precursor support exactly plane-clipped at y=.12 as an open amber back shell. No polygon selection cuts or field smoothing.','renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'size':a.size,'samples':a.samples,'device':a.device,'threads':a.threads},indent=2))
