"""Photograph real Julia thermal mineral and pressure within optical glass staging."""
import argparse, json, math, sys, hashlib, time
from pathlib import Path
import bpy, numpy as np, openvdb as vdb
from mathutils import Vector
p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path);p.add_argument('--size',type=int,default=1000);p.add_argument('--samples',type=int,default=64);p.add_argument('--step',type=int);p.add_argument('--device',default='METAL');p.add_argument('--orbit',type=float,default=.35);p.add_argument('--carrier',choices=('crystal','cutaway'),default='cutaway')
renderer_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);a.output.parent.mkdir(parents=True,exist_ok=True)
r=json.loads((a.run/'run.json').read_text());frame=r['frames'][-1] if a.step is None else next(f for f in r['frames'] if f['step']==a.step);directory=a.run/frame['directory'];shape=tuple(r['shape'])
def field(name):return np.fromfile(directory/(name+'.f32'),dtype='<f4').reshape(shape,order='F')
solid=field('solid');heat=field('temperature');pressure=field('pressure');material_pressure_max=float(np.max(np.abs(pressure)[solid>.11]))
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
s=bpy.context.scene;s.render.engine='CYCLES';s.cycles.samples=a.samples;s.cycles.use_denoising=True;s.cycles.max_bounces=16;s.cycles.transmission_bounces=12;s.cycles.volume_bounces=2
if a.device=='METAL':
 pref=bpy.context.preferences.addons['cycles'].preferences;pref.compute_device_type='METAL';pref.get_devices()
 for d in pref.devices:d.use=d.type=='METAL'
 s.cycles.device='GPU'
s.render.resolution_x=a.size;s.render.resolution_y=a.size;s.render.resolution_percentage=100;s.render.image_settings.file_format='PNG';s.world.use_nodes=True;s.world.node_tree.nodes['Background'].inputs['Color'].default_value=(.08,.12,.18,1);s.world.node_tree.nodes['Background'].inputs['Strength'].default_value=.3;s.view_settings.view_transform='AgX';s.view_settings.look='AgX - Medium High Contrast';s.view_settings.exposure=-.6

def aim(o,at=(0,0,0)):o.rotation_euler=(Vector(at)-o.location).to_track_quat('-Z','Y').to_euler()
bpy.ops.object.camera_add(location=(4.8*math.sin(a.orbit),-4.8*math.cos(a.orbit),2.1));cam=bpy.context.object;cam.data.lens=66;aim(cam,(0,0,0));s.camera=cam
for pos,power,color,size,sy in [((-3,-3,4),1150,(1,.78,.43),2,3),((2,1,2),1450,(.28,.65,1),1,3),((-.4,1,3),1700,(1,.92,.8),.4,3),((0,-2,-1),150,(1,.25,.05),2,2)]:
 bpy.ops.object.light_add(type='AREA',location=pos);o=bpy.context.object;o.data.energy=power;o.data.color=color;o.data.shape='RECTANGLE';o.data.size=size;o.data.size_y=sy;aim(o)
def mat(name,color,rough=.1,trans=0,ior=1.45,metal=0):
 m=bpy.data.materials.new(name);m.use_nodes=True;bs=m.node_tree.nodes.get('Principled BSDF')
 for k,v in [('Base Color',color),('Roughness',rough),('Transmission Weight',trans),('IOR',ior),('Metallic',metal),('Coat Weight',.25)]:bs.inputs[k].default_value=v
 return m,bs
amber,ab=mat('Pale amber carrier',(1,.85,.52,1),.035,1,1.44)
absorb=amber.node_tree.nodes.new('ShaderNodeVolumeAbsorption');absorb.inputs['Color'].default_value=(.8,.4,.09,1);absorb.inputs['Density'].default_value=.23;amber.node_tree.links.new(absorb.outputs[0],amber.node_tree.nodes.get('Material Output').inputs['Volume'])
# Clear faceted observation solid, chosen staging rather than a simulated wall.
poly=[(-.81,-.61),(-.57,-.73),(.56,-.67),(.79,-.43),(.77,.46),(.51,.64),(-.56,.68),(-.82,.39)]
if a.carrier=='cutaway':poly=[(-.82,-.65),(-.66,-.65),(-.66,.40),(.66,.40),(.66,-.65),(.82,-.65),(.82,.67),(-.82,.67)]
verts=[(x,y,z) for z in (-.79,.79) for x,y in poly];n=len(poly);faces=[tuple(reversed(range(n))),tuple(range(n,2*n))]+[(i,(i+1)%n,(i+1)%n+n,i+n) for i in range(n)]
mesh=bpy.data.meshes.new('Cut optical chamber');mesh.from_pydata(verts,[],faces);mesh.update();o=bpy.data.objects.new('Cut optical chamber',mesh);bpy.context.collection.objects.link(o);mesh.materials.append(amber);bev=o.modifiers.new('Polished edges','BEVEL');bev.width=.018;bev.segments=4;o.modifiers.new('Normals','WEIGHTED_NORMAL')
mineral,mi=mat('Copper gold mineral memory',(1,.84,.56,1),.065,1,1.20,.0)
colorattr=mineral.node_tree.nodes.new('ShaderNodeAttribute');colorattr.attribute_name='pressure_optics';mineral.node_tree.links.new(colorattr.outputs['Color'],mi.inputs['Base Color']);amplitude=mineral.node_tree.nodes.new('ShaderNodeAttribute');amplitude.attribute_name='pressure_gain';gain=mineral.node_tree.nodes.new('ShaderNodeMath');gain.operation='MULTIPLY';gain.inputs[1].default_value=3;mineral.node_tree.links.new(amplitude.outputs['Fac'],gain.inputs[0]);mi.inputs['Emission Color'].default_value=(.02,.45,1,1);mineral.node_tree.links.new(gain.outputs[0],mi.inputs['Emission Strength'])
# Two measured concentration levels: dense mineral core and a transparent dose sheath.
sheath,sh=mat('Low concentration gel interface',(1,.73,.35,1),.045,1,1.11);sh.inputs['Thin Film Thickness'].default_value=320;sh.inputs['Thin Film IOR'].default_value=1.32

def surface(data,iso,name,material):
 g=vdb.FloatGrid();g.copyFromArray(np.ascontiguousarray(data));pts,tris,quads=g.convertToPolygons(isovalue=iso,adaptivity=.025);pts=(pts+.5)*2/np.array(shape)-1
 mesh=bpy.data.meshes.new(name);mesh.from_pydata(pts.tolist(),[],tris.tolist()+quads.tolist());mesh.update();obj=bpy.data.objects.new(name,mesh);bpy.context.collection.objects.link(obj);mesh.materials.append(material)
 for f in mesh.polygons:f.use_smooth=True
 if material==mineral:
  ijk=np.clip(np.rint((pts+1)*np.array(shape)/2-.5).astype(int),0,np.array(shape)-1);sample=abs(pressure[ijk[:,0],ijk[:,1],ijk[:,2]])/max(material_pressure_max,1e-10);wave=np.clip((sample-.12)/.88,0,1);colors=np.column_stack([(1-.7*wave),(.84-.04*wave),(.56+.44*wave),np.ones(len(wave))]);attribute=mesh.color_attributes.new(name='pressure_optics',type='FLOAT_COLOR',domain='POINT');attribute.data.foreach_set('color',colors.ravel());attribute=mesh.attributes.new(name='pressure_gain',type='FLOAT',domain='POINT');attribute.data.foreach_set('value',wave.astype(np.float32))
 return {'isovalue':iso,'vertices':len(pts),'faces':len(tris)+len(quads)}
heart,heartbs=mat('Highest concentration mineral heart',(.8,.40,.07,1),.18,0,1.45,.72)
geometry=[surface(solid,.99,'Highest concentration persistent mineral',heart),surface(solid,.65,'Dense thermally converted mineral',mineral),surface(solid,.11,'Persistent low concentration sheath',sheath)]
# Saved signed pressure drives a sparse luminous volume, strongest around the
# written structure. The display gate is explicitly optical and changes no state.
pressure_scale=float(np.max(abs(pressure)));density=np.maximum(0,np.abs(pressure)/max(pressure_scale,1e-10)-.08)*np.minimum(1,solid*6)
grids=[]
for name,data in [('density',density),('temperature',heat)]:
 grid=vdb.FloatGrid();grid.name=name;grid.copyFromArray(np.ascontiguousarray(data));grid.transform=vdb.createLinearTransform(voxelSize=2/shape[0]);grid.transform.postTranslate(tuple(-1+1/shape[0] for _ in range(3)));grids.append(grid)
vpath=a.output.with_suffix('.vdb');vdb.write(str(vpath),grids=grids);data=bpy.data.volumes.new('Acoustic and thermal display');data.filepath=str(vpath.resolve());obj=bpy.data.objects.new('Acoustic and thermal display',data);bpy.context.collection.objects.link(obj)
m=bpy.data.materials.new('Sound in written mineral');m.use_nodes=True;nodes=m.node_tree.nodes;nodes.clear();out=nodes.new('ShaderNodeOutputMaterial');vol=nodes.new('ShaderNodeVolumePrincipled');vol.inputs['Color'].default_value=(.1,.55,.95,1);vol.inputs['Emission Color'].default_value=(.02,.45,1,1);vol.inputs['Density'].default_value=.2
att=nodes.new('ShaderNodeAttribute');att.attribute_name='density';gain=nodes.new('ShaderNodeMath');gain.operation='MULTIPLY';gain.inputs[1].default_value=18;m.node_tree.links.new(att.outputs['Fac'],gain.inputs[0]);m.node_tree.links.new(gain.outputs[0],vol.inputs['Emission Strength']);m.node_tree.links.new(vol.outputs['Volume'],out.inputs['Volume']);data.materials.append(m)
floor,fb=mat('Blue black photographic cove',(.018,.032,.044,1),.7);fb.inputs['Coat Weight'].default_value=0
bpy.ops.mesh.primitive_plane_add(size=200,location=(0,0,-.85));bpy.context.object.data.materials.append(floor)
s.render.filepath=str(a.output.resolve());start=time.time();bpy.ops.render.render(write_still=True)
meta={'title':r['title'],'source_run_sha256':hashlib.sha256((a.run/'run.json').read_bytes()).hexdigest(),'source_fields_sha256':{name:hashlib.sha256((directory/(name+'.f32')).read_bytes()).hexdigest() for name in ('solid','temperature','pressure')},'renderer_sha256':renderer_sha,'shape':shape,'time':frame['time'],'pressure_normalization_maximum':pressure_scale,'material_pressure_normalization_maximum':material_pressure_max,'image_sha256':hashlib.sha256(a.output.read_bytes()).hexdigest(),'geometry':geometry,'settings':{**vars(a),'run':str(a.run),'output':str(a.output)},'render_seconds':time.time()-start,'optical_mapping':'Designed faceted pale amber carrier. Native deposit isosurfaces at .99, .65 and .11 become a copper-gold heart, transparent mineral gel and thin-film transparent sheath. Absolute saved pressure / measured maximum, thresholded at .08 and gated by native mineral*6, drives cyan emission; this is a chosen display, not physical luminescence. Native pressure sampled at mineral vertices, normalized by maximum pressure within deposit>.11, also tints the core toward cyan and drives emission3× normalized amplitude above.12. No displacement or procedural interior geometry.'};a.output.with_suffix('.optical.json').write_text(json.dumps(meta,indent=2));print('DONE',meta['render_seconds'])
