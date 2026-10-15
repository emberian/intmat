"""Photograph actual deformed graph states and funded repair edges.
Skin triangles are optical panels supported only where all three mechanical
connections survive. Optical materials do not feed back into mechanics.
"""
import argparse,hashlib,json,math,sys,time
from pathlib import Path
import bpy,numpy as np
from mathutils import Vector
p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path);p.add_argument('--steps',default='all');p.add_argument('--size',type=int,default=1000);p.add_argument('--samples',type=int,default=48);p.add_argument('--device',choices=('CPU','METAL'),default='CPU');p.add_argument('--image-name',default='final.png');p.add_argument('--no-skin',action='store_true');p.add_argument('--orbit',type=float,default=.4);p.add_argument('--elevation',type=float,default=.22);p.add_argument('--distance',type=float,default=4.1);p.add_argument('--palette',choices=('pearl','cobalt'),default='pearl');p.add_argument('--target',default='0,0,0');p.add_argument('--save-scene',action='store_true')
a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);a.output.mkdir(parents=True,exist_ok=True);run=json.loads((a.run/'run.json').read_text());steps=None if a.steps=='all' else set(map(int,a.steps.split(',')));frames=[f for f in run['frames'] if steps is None or f['step'] in steps];renderer_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
scene=bpy.context.scene;scene.render.engine='CYCLES';scene.cycles.samples=a.samples;scene.cycles.use_denoising=True;scene.cycles.use_adaptive_sampling=True;scene.cycles.adaptive_threshold=.025;scene.cycles.adaptive_min_samples=min(16,a.samples);scene.cycles.max_bounces=12;scene.cycles.transmission_bounces=8
if a.device=='METAL':
 prefs=bpy.context.preferences.addons['cycles'].preferences;prefs.compute_device_type='METAL';prefs.get_devices()
 for d in prefs.devices:d.use=d.type=='METAL'
 assert any(d.use for d in prefs.devices);scene.cycles.device='GPU'
scene.render.resolution_x=a.size;scene.render.resolution_y=a.size;scene.render.resolution_percentage=100;scene.render.image_settings.file_format='PNG';scene.world.use_nodes=True;scene.world.node_tree.nodes['Background'].inputs['Color'].default_value=(.09,.065,.12,1);scene.world.node_tree.nodes['Background'].inputs['Strength'].default_value=.3;scene.view_settings.view_transform='AgX';scene.view_settings.look='AgX - Medium High Contrast';scene.view_settings.exposure=-.6

def aim(o):o.rotation_euler=(-o.location).to_track_quat('-Z','Y').to_euler()
target=Vector(tuple(map(float,a.target.split(','))));offset=Vector((a.distance*math.cos(a.elevation)*math.sin(a.orbit),-a.distance*math.cos(a.elevation)*math.cos(a.orbit),a.distance*math.sin(a.elevation)));bpy.ops.object.camera_add(location=target+offset);camera=bpy.context.object;camera.data.lens=55;camera.rotation_euler=(target-camera.location).to_track_quat('-Z','Y').to_euler();scene.camera=camera
for name,position,power,color,width,height in [('Ivory key',(-3,-3,3),900,(1,.82,.62),2,3),('Cold seam',(2,1,2),1000,(.28,.76,1),.8,3),('Rose return',(-2,1,1),750,(1,.34,.2),1.5,3),('Pearl strip',(0,-1,4),900,(.8,.9,1),.6,3)]:
 bpy.ops.object.light_add(type='AREA',location=position);o=bpy.context.object;o.name=name;o.data.energy=power;o.data.color=color;o.data.shape='RECTANGLE';o.data.size=width;o.data.size_y=height;aim(o)
def material(name,color,rough,metal=0,trans=0,ior=1.45):
 m=bpy.data.materials.new(name);m.use_nodes=True;s=m.node_tree.nodes.get('Principled BSDF')
 for key,value in [('Base Color',color),('Roughness',rough),('Metallic',metal),('Transmission Weight',trans),('IOR',ior),('Coat Weight',.3)]:s.inputs[key].default_value=value
 return m,s
old,os=material('Original oxidized woven fibres',(.12,.16,.18,1),.3,.8)
repair,rs=material('Funded copper repair fibres',(.82,.25,.055,1),.2,.78);rs.inputs['Emission Color'].default_value=(1,.12,.025,1);rs.inputs['Emission Strength'].default_value=.06
skin,ss=material('Mechanically supported pearlescent panels',(.62,.8,.83,1),.17,.15,.7,1.36)
strain=skin.node_tree.nodes.new('ShaderNodeAttribute');strain.attribute_name='strain';mult=skin.node_tree.nodes.new('ShaderNodeMath');mult.operation='MULTIPLY_ADD';mult.inputs[1].default_value=170;mult.inputs[2].default_value=260;skin.node_tree.links.new(strain.outputs['Fac'],mult.inputs[0]);skin.node_tree.links.new(mult.outputs[0],ss.inputs['Thin Film Thickness']);ss.inputs['Thin Film IOR'].default_value=1.3
if a.palette=='cobalt':
 ss.inputs['Base Color'].default_value=(.11,.31,.42,1);ss.inputs['Metallic'].default_value=.58;ss.inputs['Transmission Weight'].default_value=.16;ss.inputs['Thin Film IOR'].default_value=1.46
newskin,ns=material('New copper silk repair panels',(.94,.48,.17,1),.22,.48,.26,1.35)
floor,fs=material('Plum matte background',(.018,.012,.023,1),1,0);fs.inputs['Coat Weight'].default_value=0;fs.inputs['Specular IOR Level'].default_value=0
rows=[(-200,-1.02),(2,-1.02)]+[(2+3*math.sin(t),1.98-3*math.cos(t)) for t in np.linspace(0,math.pi/2,33)[1:]]+[(5,200)]
vertices=[(x,y,z) for y,z in rows for x in (-200,200)];faces=[(2*q,2*q+1,2*q+3,2*q+2) for q in range(len(rows)-1)];mesh=bpy.data.meshes.new('Seamless photographic cove');mesh.from_pydata(vertices,[],faces);mesh.update();back=bpy.data.objects.new('Seamless photographic cove',mesh);bpy.context.collection.objects.link(back);mesh.materials.append(floor)
for polygon in mesh.polygons:polygon.use_smooth=True

def fibres(points,edges,name,radius,mat):
 data=bpy.data.curves.new(name,'CURVE');data.dimensions='3D';data.bevel_depth=radius;data.bevel_resolution=1;data.resolution_u=1
 for i,j in edges:
  s=data.splines.new('POLY');s.points.add(1);s.points[0].co=(*points[i],1);s.points[1].co=(*points[j],1)
 o=bpy.data.objects.new(name,data);bpy.context.collection.objects.link(o);data.materials.append(mat);o['simulation_geometry']=True

def panels(points,edges,strains,generations,dimensions):
 live={tuple(sorted((int(i),int(j)))):(float(s),int(g)) for (i,j),s,g in zip(edges,strains,generations)};faces=[];repairfaces=[];vertexstrain=np.zeros(len(points))
 for (i,j),s,g in zip(edges,strains,generations):vertexstrain[i]=max(vertexstrain[i],s/.17);vertexstrain[j]=max(vertexstrain[j],s/.17)
 nx,ny,nz=dimensions
 def index(i,j,k):return i+nx*(j+ny*k)
 for k in (0,nz-1):
  for j in range(ny-1):
   for i in range(nx-1):
    v00=index(i,j,k);v10=index(i+1,j,k);v01=index(i,j+1,k);v11=index(i+1,j+1,k)
    triangles=((v00,v10,v11),(v00,v11,v01)) if (i+j)%2==0 else ((v00,v10,v01),(v10,v11,v01))
    for face in triangles:
     supporting=[live.get(tuple(sorted((face[q],face[(q+1)%3])))) for q in range(3)]
     if all(x is not None for x in supporting):faces.append(face);repairfaces.append(any(x[1]>0 for x in supporting))
 mesh=bpy.data.meshes.new('Supported optical panels');mesh.from_pydata(points.tolist(),[],faces);mesh.update();o=bpy.data.objects.new('Supported optical panels',mesh);bpy.context.collection.objects.link(o);mesh.materials.append(skin);mesh.materials.append(newskin);o['simulation_geometry']=True
 for polygon,isrepair in zip(mesh.polygons,repairfaces):polygon.use_smooth=True;polygon.material_index=int(isrepair)
 attr=mesh.attributes.new('strain','FLOAT','POINT');attr.data.foreach_set('value',np.clip(vertexstrain,0,2).astype(np.float32));return {'panels':len(faces),'repair_panels':sum(repairfaces)}
records=[]
for frame in frames:
 start=time.time();directory=a.run/frame['directory'];path=directory/'geometry.json';g=json.loads(path.read_text());count=len(g['node_ids']);points=np.asarray(g['positions']).reshape((count,3),order='F');points=points[:,[0,2,1]];points[:,1]*=-1;edges=np.asarray(g['bonds'],dtype=int).reshape((-1,2),order='F')-1;alive=np.asarray(g['damage'])==0;generation=np.asarray(g['generation']);strains=np.asarray(g['strain'])
 for o in list(bpy.data.objects):
  if o.get('simulation_geometry'):
   data=o.data;kind=o.type;bpy.data.objects.remove(o,do_unlink=True);(bpy.data.curves if kind=='CURVE' else bpy.data.meshes).remove(data)
 fibres(points,edges[alive&(generation==0)],'Original surviving fibres',.0018,old);fibres(points,edges[alive&(generation>0)],'Separately identified repair fibres',.0025,repair)
 surface=panels(points,edges[alive],strains[alive],generation[alive],g['dimensions']) if not a.no_skin else None
 output=a.output/(a.image_name if len(frames)==1 else f"frame-{frame['step']:05d}.png");scene.render.filepath=str(output.resolve())
 optical={'renderer_sha256':renderer_sha,'settings':{**vars(a),'run':str(a.run),'output':str(a.output)},'source':{'run_sha256':hashlib.sha256((a.run/'run.json').read_bytes()).hexdigest(),'geometry_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'step':frame['step'],'time':frame['time']},'camera':{'location':list(camera.location),'fixed':True,'lens_mm':55,'target':list(target)},'geometry':{'nodes':count,'live_edges':int(alive.sum()),'live_repairs':int((alive&(generation>0)).sum()),'surface':surface},'interpretation':'Exact deformed node positions, surviving original bonds and separately funded repair bonds. Thin optical panels appear only when all three mechanical connections survive. Materials, fibre radii and strain-to-film mapping are artistic; optical panels add no mechanical forces. Fixed camera; no temporal interpolation.'}
 output.with_suffix('.png.optical.json').write_text(json.dumps(optical,indent=2))
 if a.save_scene:bpy.ops.wm.save_as_mainfile(filepath=str(output.with_suffix('.blend').resolve()))
 bpy.ops.render.render(write_still=True)
 record={'step':frame['step'],'time':frame['time'],'image':output.name,'sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'render_seconds':time.time()-start};records.append(record);(a.output/'frames.json').write_text(json.dumps({'source_run':str(a.run),'fixed_camera':True,'frames':records},indent=2));print('FRAME_DONE',json.dumps(record),flush=True)
