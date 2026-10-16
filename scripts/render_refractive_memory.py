"""Photograph Julia-traced graded-index trajectories and native mineral fields.
No invented curved paths. Ray thickness/luminosity are explicit display glyphs.
"""
import argparse,hashlib,json,math,sys,time
from pathlib import Path
import bpy,numpy as np,openvdb as vdb
from mathutils import Vector
p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('optics',type=Path);p.add_argument('output',type=Path)
p.add_argument('--size',type=int,default=640);p.add_argument('--samples',type=int,default=32);p.add_argument('--device',choices=('CPU','METAL'),default='METAL');p.add_argument('--radius',type=float,default=.0010);p.add_argument('--gain',type=float,default=2.5);p.add_argument('--orbit',type=float,default=.50);p.add_argument('--elevation',type=float,default=.20);p.add_argument('--distance',type=float,default=5.6);p.add_argument('--material',choices=('volume','surface','none'),default='volume');p.add_argument('--screen',action='store_true');p.add_argument('--save-scene',action='store_true');p.add_argument('--wavelength',type=int,default=0);p.add_argument('--ray-stride',type=int,default=1);p.add_argument('--display-zmin',type=float,default=-1.08);p.add_argument('--display-zmax',type=float,default=1.45);p.add_argument('--no-plate',action='store_true');p.add_argument('--launch-radius',type=float,default=10.);p.add_argument('--ray-palette',choices=('spectral','gold','ice'),default='spectral');p.add_argument('--no-floor',action='store_true');p.add_argument('--surface-style',choices=('ghost','glass'),default='ghost');p.add_argument('--skin-iso',type=float,default=.16);p.add_argument('--core-iso',type=float,default=0);p.add_argument('--cut-y',type=float,default=-2.);p.add_argument('--rays-in-material',type=float,default=0);p.add_argument('--skin-absorption',type=float,default=1.2);p.add_argument('--slate',action='store_true')
a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);a.output.parent.mkdir(parents=True,exist_ok=True)
run=json.loads((a.run/'run.json').read_text());optical=json.loads((a.optics/'observation.json').read_text());frame=next(f for f in run['frames'] if f['step']==optical['step']);shape=tuple(run['shape'])
solid=np.fromfile(a.run/frame['directory']/'solid.f32',dtype='<f4').reshape(shape,order='F')
assert hashlib.sha256((a.run/frame['directory']/'solid.f32').read_bytes()).hexdigest()==optical['source_solid_sha256']
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
scene=bpy.context.scene;scene.render.engine='CYCLES';scene.cycles.samples=a.samples;scene.cycles.use_denoising=True;scene.cycles.use_adaptive_sampling=True;scene.cycles.adaptive_threshold=.025;scene.cycles.max_bounces=10;scene.cycles.volume_bounces=2
if a.device=='METAL':
 prefs=bpy.context.preferences.addons['cycles'].preferences;prefs.compute_device_type='METAL';prefs.get_devices()
 for d in prefs.devices:d.use=d.type=='METAL'
 assert any(d.use for d in prefs.devices);scene.cycles.device='GPU'
else:scene.cycles.device='CPU'
scene.render.resolution_x=a.size;scene.render.resolution_y=a.size;scene.render.resolution_percentage=100;scene.render.image_settings.file_format='PNG';scene.render.filepath=str(a.output.resolve())
scene.world.use_nodes=True;scene.world.node_tree.nodes['Background'].inputs['Color'].default_value=(.007,.011,.022,1);scene.world.node_tree.nodes['Background'].inputs['Strength'].default_value=.16
if a.slate:
 scene.world.node_tree.nodes['Background'].inputs['Color'].default_value=(.055,.085,.13,1);scene.world.node_tree.nodes['Background'].inputs['Strength'].default_value=.5
scene.view_settings.view_transform='AgX';scene.view_settings.look='AgX - Medium High Contrast';scene.view_settings.exposure=-.2

def aim(o,target=(0,0,.15)):o.rotation_euler=(Vector(target)-o.location).to_track_quat('-Z','Y').to_euler()
bpy.ops.object.camera_add(location=(a.distance*math.cos(a.elevation)*math.sin(a.orbit),-a.distance*math.cos(a.elevation)*math.cos(a.orbit),a.distance*math.sin(a.elevation)+.15));camera=bpy.context.object;camera.data.lens=57;aim(camera);scene.camera=camera
for pos,power,color,sx,sy in [((-3,-2,3),650,(.45,.75,1),1,4),((2,1,3),750,(1,.55,.23),1,3),((-2,2,.2),400,(.65,.9,1),1,2)]:
 bpy.ops.object.light_add(type='AREA',location=pos);o=bpy.context.object;o.data.energy=power;o.data.color=color;o.data.shape='RECTANGLE';o.data.size=sx;o.data.size_y=sy;aim(o)

def mat(name,color,rough=.3,metal=0):
 m=bpy.data.materials.new(name);m.use_nodes=True;s=m.node_tree.nodes.get('Principled BSDF');s.inputs['Base Color'].default_value=(*color,1);s.inputs['Roughness'].default_value=rough;s.inputs['Metallic'].default_value=metal;return m,s
black,_=mat('Midnight optical bench',(.004,.007,.012),.36,.20)
bpy.ops.mesh.primitive_plane_add(size=200,location=(0,0,-1.13));bpy.context.object.data.materials.append(black)
if a.no_floor:bpy.data.objects.remove(bpy.context.object,do_unlink=True)
# A thin launch plate is optical staging, outside the simulated material box.
bpy.ops.mesh.primitive_cube_add(size=1,location=(0,0,-1.105));plate=bpy.context.object;plate.name='Launch aperture / display apparatus';plate.scale=(1.28,1.28,.035);plate.data.materials.append(black)
bevel=plate.modifiers.new('Machined edge','BEVEL');bevel.width=.02;bevel.segments=3
if a.no_plate:bpy.data.objects.remove(plate,do_unlink=True)
colors={460:(.055,.26,1.),550:(.11,1.,.43),630:(1.,.22,.035)}
def sample_material(points):
 pos=np.clip((points+1)*np.array(shape)/2-.5,0,np.array(shape)-1);lo=np.floor(pos).astype(int);hi=np.minimum(lo+1,np.array(shape)-1);frac=pos-lo;value=np.zeros(len(points),dtype=np.float32)
 for x in (0,1):
  for y in (0,1):
   for z in (0,1):
    bits=(x,y,z);idx=[hi[:,d] if bits[d] else lo[:,d] for d in range(3)];weight=np.prod([frac[:,d] if bits[d] else 1-frac[:,d] for d in range(3)],axis=0);value+=solid[tuple(idx)]*weight
 value[np.any(np.abs(points)>1,axis=1)]=0
 return value
counts={}
for wavelength,color in colors.items():
 if a.wavelength and wavelength!=a.wavelength:continue
 color=color if a.ray_palette=='spectral' else (1.,.40,.08) if a.ray_palette=='gold' else (.10,.62,1.)
 material,bs=mat(f'{wavelength} nm ray display',color,.24,.25);bs.inputs['Emission Color'].default_value=(*color,1);bs.inputs['Emission Strength'].default_value=a.gain
 curve=bpy.data.curves.new(f'Recorded {wavelength} nm paths','CURVE');curve.dimensions='3D';curve.resolution_u=1;curve.bevel_depth=a.radius;curve.bevel_resolution=2
 rays=[r for r in optical['rays'] if r['wavelength']==wavelength];rays=[r for i,r in enumerate(rays) if i%optical['visual_grid']%a.ray_stride==0 and (i//optical['visual_grid'])%a.ray_stride==0];rays=[r for r in rays if math.hypot(*r['launch'][:2])<=a.launch_radius];counts[wavelength]=len(rays)
 for r in rays:
  full=np.asarray(r['path'],dtype=np.float32);kept=np.flatnonzero((full[:,2]>=a.display_zmin)&(full[:,2]<=a.display_zmax)&(sample_material(full)>=a.rays_in_material))
  for indices in np.split(kept,np.flatnonzero(np.diff(kept)>1)+1):
   if len(indices)<2:continue
   xyz=full[indices];spline=curve.splines.new('POLY');spline.points.add(len(xyz)-1);spline.points.foreach_set('co',np.column_stack([xyz,np.ones(len(xyz))]).ravel())
 obj=bpy.data.objects.new(f'Traced {wavelength} nm paths',curve);bpy.context.collection.objects.link(obj);curve.materials.append(material)
if a.material=='volume':
 grid=vdb.FloatGrid();grid.name='solid';grid.copyFromArray(np.ascontiguousarray(solid));grid.transform=vdb.createLinearTransform(voxelSize=2/shape[0]);grid.transform.postTranslate((-1+1/shape[0],)*3)
 path=a.output.with_suffix('.vdb');vdb.write(str(path),grids=[grid]);volume=bpy.data.volumes.new('Actual mineral concentration');volume.filepath=str(path.resolve());obj=bpy.data.objects.new('Finite mineral / optical haze',volume);bpy.context.collection.objects.link(obj)
 m=bpy.data.materials.new('Mineral observation haze');m.use_nodes=True;nodes=m.node_tree.nodes;nodes.clear();out=nodes.new('ShaderNodeOutputMaterial');attribute=nodes.new('ShaderNodeAttribute');attribute.attribute_name='solid';density=nodes.new('ShaderNodeMath');density.operation='MULTIPLY';density.inputs[1].default_value=.7;m.node_tree.links.new(attribute.outputs['Fac'],density.inputs[0]);scatter=nodes.new('ShaderNodeVolumePrincipled');scatter.inputs['Color'].default_value=(.3,.6,.68,1);m.node_tree.links.new(density.outputs[0],scatter.inputs['Density']);m.node_tree.links.new(scatter.outputs[0],out.inputs['Volume']);volume.materials.append(m)
elif a.material=='surface':
 observed=solid.copy();ys=(np.arange(shape[1])+.5)*2/shape[1]-1;observed[:,ys<a.cut_y,:]=0
 def contour(iso,name,material):
  grid=vdb.FloatGrid();grid.copyFromArray(np.ascontiguousarray(observed));points,tris,quads=grid.convertToPolygons(isovalue=iso,adaptivity=.02);points=(points+.5)*2/np.array(shape)-1
  mesh=bpy.data.meshes.new(name);mesh.from_pydata(points.tolist(),[],tris.tolist()+quads.tolist());mesh.update();obj=bpy.data.objects.new(name,mesh);bpy.context.collection.objects.link(obj);mesh.materials.append(material)
  for face in mesh.polygons:face.use_smooth=True
  return obj
 if a.surface_style=='ghost':
  m,bs=mat('Mineral ghost display',(.15,.38,.43),.22,.3);transparent=m.node_tree.nodes.new('ShaderNodeBsdfTransparent');mix=m.node_tree.nodes.new('ShaderNodeMixShader');mix.inputs[0].default_value=.88;m.node_tree.links.new(bs.outputs[0],mix.inputs[1]);m.node_tree.links.new(transparent.outputs[0],mix.inputs[2]);m.node_tree.links.new(mix.outputs[0],m.node_tree.nodes.get('Material Output').inputs[0])
 else:
  m,bs=mat('Warm refractive mineral contour',(1,.78,.46),.075);bs.inputs['Transmission Weight'].default_value=.98;bs.inputs['IOR'].default_value=1.34;bs.inputs['Coat Weight'].default_value=.25
  absorption=m.node_tree.nodes.new('ShaderNodeVolumeAbsorption');absorption.inputs['Color'].default_value=(.65,.22,.045,1);absorption.inputs['Density'].default_value=a.skin_absorption;m.node_tree.links.new(absorption.outputs[0],m.node_tree.nodes.get('Material Output').inputs['Volume'])
 contour(a.skin_iso,f'Saved solid {a.skin_iso} optical boundary',m)
 if a.core_iso>0:
  m,bs=mat('Dense mineral / plum gold inclusion',(.36,.13,.065),.28,.7);bs.inputs['Coat Weight'].default_value=.2
  contour(a.core_iso,f'Saved solid {a.core_iso} dense contour',m)
# The count screen is retained as a separate actual-observation image even when
# omitted from the sculpture composition. Display blur is a fixed pixel kernel.
N=optical['screen_pixels'];rgb=np.zeros((N,N,3),dtype=np.float32)
for screen in optical['screens']:
 h=np.array(screen['counts'],dtype=np.float32).reshape((N,N),order='F')*N*N/screen['total']
 for axis in (0,1):
  pad=np.pad(h,[(2,2) if d==axis else (0,0) for d in (0,1)])
  h=sum(weight*np.take(pad,range(shift,shift+N),axis=axis) for shift,weight in enumerate([1/16,4/16,6/16,4/16,1/16]))
 rgb+=(1-np.exp(-.16*h))[...,None]*np.array(colors[int(screen['wavelength'])])/2
pixels=np.concatenate([np.transpose(rgb,(1,0,2)),np.ones((N,N,1),dtype=np.float32)],axis=2)
image=bpy.data.images.new('Equal-weight traced flux',width=N,height=N,float_buffer=True);image.pixels.foreach_set(pixels.ravel());image.filepath_raw=str(a.output.with_suffix('.screen.png').resolve());image.file_format='PNG';image.save()
if a.screen:
 material=bpy.data.materials.new('Receiving plane / computed ray density');material.use_nodes=True;nodes=material.node_tree.nodes;nodes.clear();out=nodes.new('ShaderNodeOutputMaterial');emission=nodes.new('ShaderNodeEmission');tex=nodes.new('ShaderNodeTexImage');tex.image=image;material.node_tree.links.new(tex.outputs['Color'],emission.inputs['Color']);emission.inputs['Strength'].default_value=1.4;material.node_tree.links.new(emission.outputs[0],out.inputs['Surface'])
 bpy.ops.mesh.primitive_plane_add(size=optical['screen_span'],location=(0,0,optical['target_z']));bpy.context.object.name='Actual target z plane';bpy.context.object.data.materials.append(material)
metadata=dict(renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),observation_sha256=hashlib.sha256((a.optics/'observation.json').read_bytes()).hexdigest(),settings={**vars(a),'run':str(a.run),'optics':str(a.optics),'output':str(a.output)},source=dict(shape=shape,time=optical['time'],step=optical['step'],checks=optical['checks']),ray_counts=counts,interpretation='Julia GRIN path polylines, equal-weight input rays. Luminous tubes are visible path glyphs, not computed scattering radiance. Mineral haze/ghost contour and dark launch plate are optical staging. Count screen uses fixed five-tap display blur and1-exp(-.16 normalized count); no interference/Fresnel or coherent-wave caustic calculation.')
a.output.with_suffix('.png.optical.json').write_text(json.dumps(metadata,indent=2))
if a.save_scene:bpy.ops.wm.save_as_mainfile(filepath=str(a.output.with_suffix('.blend').resolve()))
bpy.ops.render.render(write_still=True)
print('RENDER_DONE',a.output,flush=True)
