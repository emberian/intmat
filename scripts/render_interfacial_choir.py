"""Photograph saved Julia Nacre Fault fields as interacting optical phases.
No shape key, procedural displacement, time interpolation or camera animation.
Blender: -b --python this.py -- RUN OUT --steps 0,20,...
"""
import argparse, hashlib, json, math, os, sys, time
from pathlib import Path
import bpy, numpy as np, openvdb as vdb
from mathutils import Vector
p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path)
p.add_argument('--device',choices=('CPU','METAL'),default='METAL');p.add_argument('--image-name',default='final.png');p.add_argument('--steps',default='all');p.add_argument('--size',type=int,default=720);p.add_argument('--samples',type=int,default=32)
p.add_argument('--solid-iso',type=float,default=.16);p.add_argument('--core-iso',type=float,default=.60);p.add_argument('--pressure-iso',type=float,default=.12)
p.add_argument('--skin-transmission',type=float,default=.9);p.add_argument('--skin-roughness',type=float,default=.14)
p.add_argument('--no-wave',action='store_true');p.add_argument('--save-scene',action='store_true');p.add_argument('--orbit',type=float,default=.35);p.add_argument('--elevation',type=float,default=.32);p.add_argument('--distance',type=float,default=4.0)
renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);a.output.mkdir(parents=True,exist_ok=True)
run=json.loads((a.run/'run.json').read_text());frames=run['frames'];steps=None if a.steps=='all' else set(map(int,a.steps.split(',')));frames=[f for f in frames if steps is None or f['step'] in steps]
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
scene=bpy.context.scene;scene.render.engine='CYCLES';scene.cycles.samples=a.samples;scene.cycles.use_denoising=True
scene.cycles.use_adaptive_sampling=True;scene.cycles.adaptive_threshold=.025;scene.cycles.adaptive_min_samples=min(16,a.samples)
scene.cycles.max_bounces=14;scene.cycles.transmission_bounces=10;scene.cycles.transparent_max_bounces=12
if a.device=='METAL':
    prefs=bpy.context.preferences.addons['cycles'].preferences;prefs.compute_device_type='METAL';prefs.get_devices()
    for device in prefs.devices:device.use=device.type=='METAL'
    assert any(d.use for d in prefs.devices);scene.cycles.device='GPU'
else:scene.cycles.device='CPU'
scene.render.resolution_x=a.size;scene.render.resolution_y=a.size;scene.render.resolution_percentage=100
scene.render.image_settings.file_format='PNG';scene.world.use_nodes=True
scene.world.node_tree.nodes['Background'].inputs['Color'].default_value=(.045,.065,.09,1)
scene.world.node_tree.nodes['Background'].inputs['Strength'].default_value=.25
scene.view_settings.view_transform='AgX';scene.view_settings.look='AgX - Medium High Contrast';scene.view_settings.exposure=-.65

def aim(o,target=(0,0,0)):o.rotation_euler=(Vector(target)-o.location).to_track_quat('-Z','Y').to_euler()
bpy.ops.object.camera_add(location=(a.distance*math.cos(a.elevation)*math.sin(a.orbit),-a.distance*math.cos(a.elevation)*math.cos(a.orbit),a.distance*math.sin(a.elevation)))
camera=bpy.context.object;camera.data.lens=55;aim(camera);scene.camera=camera
for name,pos,power,color,size,size_y in [('Ivory ribbon',(-3,-4,4),650,(1,.9,.72),1.2,3.8),('Petrol edge',(3,1,2),800,(.3,.72,1),1.2,3.2),('Warm seam',(-2,1,1.2),420,(1,.26,.08),2,3),('Top slit',(0,0,4),700,(.72,.85,1),.7,3)]:
    bpy.ops.object.light_add(type='AREA',location=pos);light=bpy.context.object;light.name=name;light.data.energy=power;light.data.color=color;light.data.shape='RECTANGLE';light.data.size=size;light.data.size_y=size_y;aim(light)
bpy.ops.mesh.primitive_plane_add(size=200,location=(0,0,-1.03));ground=bpy.context.object;ground.name='Studio ground'
m=bpy.data.materials.new('Midnight satin');m.use_nodes=True;s=m.node_tree.nodes.get('Principled BSDF');s.inputs['Base Color'].default_value=(.004,.007,.012,1);s.inputs['Roughness'].default_value=.55;s.inputs['Metallic'].default_value=.08;ground.data.materials.append(m)

def material(name,color,rough,metal=0,transmission=0,ior=1.46):
    m=bpy.data.materials.new(name);m.use_nodes=True;s=m.node_tree.nodes.get('Principled BSDF')
    for key,value in [('Base Color',color),('Roughness',rough),('Metallic',metal),('Transmission Weight',transmission),('IOR',ior),('Coat Weight',.3),('Coat Roughness',.06)]:s.inputs[key].default_value=value
    return m,s
shell,ss=material('Pressure-tuned nacre',(.25,.65,.63,1),a.skin_roughness,0,a.skin_transmission,1.34)
# Film thickness and color depend on saved fields sampled at the actual surface.
att=shell.node_tree.nodes.new('ShaderNodeAttribute');att.attribute_name='material_activity'
ramp=shell.node_tree.nodes.new('ShaderNodeValToRGB');ramp.color_ramp.elements.remove(ramp.color_ramp.elements[1])
for ix,(pos,col) in enumerate([(0,(.45,.72,.83,1)),(.35,(.62,.91,.84,1)),(.62,(.88,.68,.85,1)),(1,(.98,.88,.64,1))]):
    e=ramp.color_ramp.elements[0] if ix==0 else ramp.color_ramp.elements.new(pos);e.position=pos;e.color=col
shell.node_tree.links.new(att.outputs['Fac'],ramp.inputs['Fac']);shell.node_tree.links.new(ramp.outputs['Color'],ss.inputs['Base Color'])
film=shell.node_tree.nodes.new('ShaderNodeAttribute');film.attribute_name='film_nm';shell.node_tree.links.new(film.outputs['Fac'],ss.inputs['Thin Film Thickness']);ss.inputs['Thin Film IOR'].default_value=1.34
core,cs=material('Dense mineral / warm palladium',(.63,.25,.065,1),.34,.78,.0)
cs.inputs['Anisotropic'].default_value=.35
wave,ws=material('Pressure crest optical surrogate',(.65,.86,1,1),.08,.0,.94,1.12)
ws.inputs['Thin Film Thickness'].default_value=390;ws.inputs['Thin Film IOR'].default_value=1.28
ws.inputs['Emission Color'].default_value=(.1,.35,.48,1);ws.inputs['Emission Strength'].default_value=.035
# Fade the observation boundary before its discretized face cut. This is an
# explicitly optical visibility transfer; saved pressure and geometry are intact.
visibility=wave.node_tree.nodes.new('ShaderNodeAttribute');visibility.attribute_name='crop_opacity'
transparent=wave.node_tree.nodes.new('ShaderNodeBsdfTransparent')
blend=wave.node_tree.nodes.new('ShaderNodeMixShader')
wave.node_tree.links.new(visibility.outputs['Fac'],blend.inputs[0])
wave.node_tree.links.new(transparent.outputs[0],blend.inputs[1])
wave.node_tree.links.new(ws.outputs[0],blend.inputs[2])
wave.node_tree.links.new(blend.outputs[0],wave.node_tree.nodes.get('Material Output').inputs['Surface'])

def sample(field,points):
    n=np.array(field.shape);pos=np.clip((points+1)*n/2-.5,0,n-1);lo=np.floor(pos).astype(int);hi=np.minimum(lo+1,n-1);f=pos-lo;result=np.zeros(len(points),dtype=np.float32)
    for x in (0,1):
      for y in (0,1):
       for z in (0,1):
        indices=[hi[:,d] if bit else lo[:,d] for d,bit in enumerate((x,y,z))]
        weight=np.prod([f[:,d] if bit else 1-f[:,d] for d,bit in enumerate((x,y,z))],axis=0)
        result+=field[tuple(indices)]*weight
    return result

def surface(field,iso,name,mat,attributes=None,clip=None):
    grid=vdb.FloatGrid();grid.copyFromArray(np.ascontiguousarray(field));points,tris,quads=grid.convertToPolygons(isovalue=iso,adaptivity=.035)
    points=(points+.5)*2/np.array(field.shape)-1
    faces=tris.tolist()+quads.tolist()
    if clip is not None:
        keep=clip(points);faces=[f for f in faces if all(keep[i] for i in f)]
    mesh=bpy.data.meshes.new(name);mesh.from_pydata(points.tolist(),[],faces);mesh.update();o=bpy.data.objects.new(name,mesh);bpy.context.collection.objects.link(o);o.data.materials.append(mat)
    for poly in mesh.polygons:poly.use_smooth=True
    if attributes:
        for key,values in attributes.items():
            attr=mesh.attributes.new(key,'FLOAT','POINT');attr.data.foreach_set('value',values(points).astype(np.float32))
    return o,{'vertices':len(points),'faces':len(faces),'isovalue':iso,'field':name}

records=[]
for frame in frames:
    start=time.time();directory=a.run/frame['directory'];meta=json.loads((directory/'volume.json').read_text());shape=tuple(meta['shape'])
    fields={name:np.fromfile(directory/(name+'.f32'),dtype='<f4').reshape(shape,order='F') for name in ('solid','activity','pressure')}
    for obj in list(bpy.data.objects):
        if obj.get('simulation_geometry'):mesh=obj.data;bpy.data.objects.remove(obj,do_unlink=True);bpy.data.meshes.remove(mesh)
    shellobj,shellmeta=surface(fields['solid'],a.solid_iso,'solid / translucent skin',shell,{'material_activity':lambda points:np.clip((sample(fields['activity'],points)+1)*.5,0,1),'film_nm':lambda points:380+170*np.tanh(4*sample(fields['pressure'],points))})
    coreobj,coremeta=surface(fields['solid'],a.core_iso,'solid / dense core',core)
    objects=[shellobj,coreobj];meshes=[shellmeta,coremeta]
    if not a.no_wave:
        o,m=surface(fields['pressure'],a.pressure_iso,'pressure / positive crest',wave,attributes={'crop_opacity':lambda points:np.minimum(np.clip((.86-np.linalg.norm(points[:,[0,2]],axis=1))/.07,0,1),np.clip((.58-np.abs(points[:,1]))/.05,0,1))**3},clip=lambda points:(np.linalg.norm(points[:,[0,2]],axis=1)<.86)&(np.abs(points[:,1])<.58));objects.append(o);meshes.append(m)
    for obj in objects:obj['simulation_geometry']=True
    output=a.output/(a.image_name if len(frames)==1 else f"frame-{frame['step']:05d}.png")
    scene.render.filepath=str(output.resolve())
    optical={'renderer_sha256':renderer_sha256,'settings':{**vars(a),'run':str(a.run),'output':str(a.output)},'source':{'shape':shape,'step':frame['step'],'time':frame['time'],'run_sha256':hashlib.sha256((a.run/'run.json').read_bytes()).hexdigest()},'meshes':meshes,'camera':{'location':list(camera.location),'lens_mm':55,'fixed':True},'interpretation':{'skin':'Saved solid isosurface; artistic translucent nacre; color from saved recurrent activity; thin-film thickness 380+170*tanh(4*pressure) nm is optical mapping, not simulated film thickness.','core':'Higher solid isosurface with chosen warm metallic BSDF.','wave':'Positive saved pressure isosurface with radius .86 / y .58 display clipping, cubic opacity fade over the outer .07 radius / .05 y widths; refractive optical surrogate, not a material phase.','lighting':'Fixed area lights; no generated texture or displacement.'}}
    output.with_suffix('.png.optical.json').write_text(json.dumps(optical,indent=2))
    if a.save_scene:bpy.ops.wm.save_as_mainfile(filepath=str(output.with_suffix('.blend').resolve()))
    bpy.ops.render.render(write_still=True)
    record={'step':frame['step'],'time':frame['time'],'image':output.name,'sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'render_seconds':time.time()-start,'meshes':meshes};records.append(record)
    (a.output/'frames.json').write_text(json.dumps({'source_run':str(a.run),'shape':shape,'fixed_camera':True,'frames':records},indent=2));print('FRAME_DONE',json.dumps(record),flush=True)
