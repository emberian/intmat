"""Photograph computed mineral/current/damage inside a designed colored optical prism."""
import argparse,json,sys,hashlib,math
from pathlib import Path
import bpy,numpy as np,openvdb as vdb
from mathutils import Vector
p=argparse.ArgumentParser();p.add_argument('frame',type=Path);p.add_argument('output',type=Path);p.add_argument('--size',type=int,default=640);p.add_argument('--samples',type=int,default=16);p.add_argument('--iso',type=float,default=.22);p.add_argument('--device',choices=('CPU','METAL'),default='CPU');p.add_argument('--threads',type=int,default=2);a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
if a.threads<1:p.error('--threads must be positive')
meta=json.loads((a.frame/'volume.json').read_text());shape=tuple(meta['shape']);nodal='sample_convention' in meta;h=2/(shape[0]-1) if nodal else 2/shape[0]
fields={name:np.fromfile(a.frame/(name+'.f32'),dtype='<f4').reshape(shape,order='F') for name in ['deposit','activity','damage']}
grid=vdb.FloatGrid();grid.copyFromArray(np.ascontiguousarray(fields['deposit']));points,tri,quad=grid.convertToPolygons(isovalue=a.iso,adaptivity=.05)
indices=np.clip(np.rint(points).astype(int),0,np.array(shape)-1);points=points*h-1 if nodal else (points+.5)*h-1
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
mesh=bpy.data.meshes.new('Computed mineral level set');mesh.from_pydata(points.tolist(),[],tri.tolist()+quad.tolist());mesh.update();obj=bpy.data.objects.new('Conductive mineral / actual3D state',mesh);bpy.context.collection.objects.link(obj)
for poly in mesh.polygons:poly.use_smooth=True
for name,field in [('joule',fields['activity']),('damage',fields['damage'])]:
 attr=mesh.attributes.new(name,'FLOAT','POINT');attr.data.foreach_set('value',field[tuple(indices.T)].astype(np.float32))
def mat(name,color,rough,trans,ior):
 m=bpy.data.materials.new(name);m.use_nodes=True;s=m.node_tree.nodes.get('Principled BSDF');s.inputs['Base Color'].default_value=color;s.inputs['Roughness'].default_value=rough;s.inputs['Transmission Weight'].default_value=trans;s.inputs['IOR'].default_value=ior;s.inputs['Coat Weight'].default_value=.35;return m,s
m,s=mat('Opaline conductor / current and damage optics',(.57,.86,.75,1),.13,.47,1.39);s.inputs['Metallic'].default_value=.12;s.inputs['Thin Film IOR'].default_value=1.33;s.inputs['Emission Color'].default_value=(1,.045,.22,1)
nodes=m.node_tree.nodes;links=m.node_tree.links
for name,target,scale,offset in [('joule','Emission Strength',5.,0),('damage','Roughness',.35,.08),('damage','Thin Film Thickness',450.,250.)]:
 attr=nodes.new('ShaderNodeAttribute');attr.attribute_name=name;gain=nodes.new('ShaderNodeMath');gain.operation='MULTIPLY_ADD';gain.inputs[1].default_value=scale;gain.inputs[2].default_value=offset;links.new(attr.outputs['Fac'],gain.inputs[0]);links.new(gain.outputs[0],s.inputs[target])
mesh.materials.append(m)
# The prism is photographic staging; all of its internal mineral comes from the state.
medium,ms=mat('Blue optical medium / designed carrier',(.49,.75,.94,1),.035,1,1.33);absorb=medium.node_tree.nodes.new('ShaderNodeVolumeAbsorption');absorb.inputs['Color'].default_value=(.14,.43,.68,1);absorb.inputs['Density'].default_value=.22;medium.node_tree.links.new(absorb.outputs[0],medium.node_tree.nodes.get('Material Output').inputs['Volume'])
poly=[(-1.09,-.86),(-.88,-1.07),(.82,-1.07),(1.09,-.80),(1.09,.84),(.86,1.07),(-.85,1.07),(-1.09,.83)]
verts=[(x,y,z) for y in (-.94,.94) for x,z in poly];count=len(poly);faces=[tuple(reversed(range(count))),tuple(range(count,2*count))]+[(i,(i+1)%count,(i+1)%count+count,i+count) for i in range(count)]
faces=[tuple(reversed(f)) for f in faces]
me=bpy.data.meshes.new('Designed optical prism');me.from_pydata(verts,[],faces);me.update();o=bpy.data.objects.new('Colored observation medium',me);bpy.context.collection.objects.link(o);me.materials.append(medium);bevel=o.modifiers.new('Polished optical edge','BEVEL');bevel.width=.018;bevel.segments=3;o.modifiers.new('Weighted normals','WEIGHTED_NORMAL')
sc=bpy.context.scene;sc.render.engine='CYCLES';sc.cycles.device='CPU';sc.render.threads_mode='FIXED';sc.render.threads=a.threads;
if a.device=='METAL':
 prefs=bpy.context.preferences.addons['cycles'].preferences;prefs.compute_device_type='METAL';prefs.get_devices()
 for device in prefs.devices:device.use=device.type=='METAL'
 if not any(device.use for device in prefs.devices):raise RuntimeError('No enabled Metal device')
 sc.cycles.device='GPU'
sc.cycles.samples=a.samples;sc.cycles.use_denoising=True;sc.cycles.max_bounces=16;sc.cycles.transmission_bounces=12
sc.render.resolution_x=a.size;sc.render.resolution_y=a.size;sc.render.resolution_percentage=100;sc.world.use_nodes=True;sc.world.node_tree.nodes['Background'].inputs[0].default_value=(.10,.16,.23,1);sc.world.node_tree.nodes['Background'].inputs[1].default_value=.3;sc.view_settings.view_transform='AgX';sc.view_settings.look='AgX - Medium High Contrast';sc.view_settings.exposure=.1
def aim(o,at=(0,0,0)):o.rotation_euler=(Vector(at)-o.location).to_track_quat('-Z','Y').to_euler()
bpy.ops.object.camera_add(location=(3,-4,2.6));camera=bpy.context.object;camera.data.type='ORTHO';camera.data.ortho_scale=3.65;aim(camera);sc.camera=camera
for name,pos,power,color,size,sy in [('Warm silk',(-3,-3,4),750,(1,.76,.51),3,3),('White edge',(3,1,3),1100,(.76,.9,1),.5,3),('Rose ribbon',(-2,1,.2),450,(1,.25,.45),.6,2.8)]:
 bpy.ops.object.light_add(type='AREA',location=pos);o=bpy.context.object;o.name=name;o.data.energy=power;o.data.color=color;o.data.shape='RECTANGLE';o.data.size=size;o.data.size_y=sy;aim(o)
floor,fs=mat('Slate photographic cove',(.025,.044,.061,1),.65,0,1.45);bpy.ops.mesh.primitive_plane_add(size=200,location=(0,0,-1.1));bpy.context.object.data.materials.append(floor)
a.output.parent.mkdir(parents=True,exist_ok=True);sc.render.image_settings.file_format='PNG';sc.render.filepath=str(a.output.resolve());bpy.ops.render.render(write_still=True)
a.output.with_suffix('.optical.json').write_text(json.dumps({'source_frame':a.frame.name,'fields':'Actual deposit level set; actual local Joule power drives pink emission; actual damage drives roughness and thin-film optics.','source_shape':shape,'node_centered':nodal,'isovalue':a.iso,'vertices':len(points),'polygons':len(tri)+len(quad),'carrier':'Designed blue absorbing polished prism; optical staging, not simulated material geometry.','renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'size':a.size,'samples':a.samples,'device':a.device,'threads':a.threads},indent=2))
