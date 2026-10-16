"""Observe actual voltage iso-surfaces through surviving mechanical tetrahedra."""
import bpy,bmesh,numpy as np,json,argparse,sys,hashlib,itertools
from pathlib import Path
from mathutils import Vector
p=argparse.ArgumentParser();p.add_argument('state',type=Path);p.add_argument('output',type=Path);p.add_argument('--size',type=int,default=640);p.add_argument('--samples',type=int,default=16);p.add_argument('--device',choices=['CPU','METAL'],default='CPU');a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
s=json.loads(a.state.read_text());nx,ny,nz=s['dims']
def columns(key,width):
 v=np.array(s[key]);return v.reshape(-1,width) if v.ndim==1 else v
pos=columns('position',3);edges=columns('bonds',2).astype(int)-1;damage=np.array(s['damage']);voltage=np.array(s['potential']);dose=np.array(s['dose']);edge={tuple(sorted(e)):i for i,e in enumerate(edges)}
tetpattern=((0,1,3,7),(0,3,2,7),(0,2,6,7),(0,6,4,7),(0,4,5,7),(0,5,1,7));tets=[];boundary={};removed=0
for k in range(nz-1):
 for j in range(ny-1):
  for i in range(nx-1):
   q=i+nx*j+nx*ny*k;cube=[q,q+1,q+nx,q+nx+1,q+nx*ny,q+nx*ny+1,q+nx*ny+nx,q+nx*ny+nx+1]
   for pattern in tetpattern:
    ids=[cube[v] for v in pattern]
    if any(damage[edge[tuple(sorted(pair))]]>.5 for pair in itertools.combinations(ids,2)):removed+=1;continue
    tets.append(ids)
    for f in ((ids[0],ids[2],ids[1]),(ids[0],ids[1],ids[3]),(ids[0],ids[3],ids[2]),(ids[1],ids[2],ids[3])):
     key=tuple(sorted(f));boundary[key]=None if key in boundary else f
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
def material(name,color,rough,metal,trans,ior):
 m=bpy.data.materials.new(name);m.use_nodes=True;q=m.node_tree.nodes.get('Principled BSDF');q.inputs['Base Color'].default_value=color;q.inputs['Roughness'].default_value=rough;q.inputs['Metallic'].default_value=metal;q.inputs['Transmission Weight'].default_value=trans;q.inputs['IOR'].default_value=ior;q.inputs['Coat Weight'].default_value=.18;return m,q
mats=[]
for j in range(16):
 t=j/15;m,q=material('Stored electrical dose '+str(j),(.13+.39*t,.085+.19*t,.16-.11*t,1),.35-.1*t,.55,.15,1.45);mats.append(m)
vertices=[];faces=[];materials=[];cache={};isos=[.12,.22,.32,.42,.52,.62,.72,.82,.90]
for level in isos:
 for ids in tets:
  if not min(voltage[ids])<level<max(voltage[ids]):continue
  polygon=[];direction=None
  for i,j in itertools.combinations(ids,2):
   if (voltage[i]<level)!=(voltage[j]<level):
    if direction is None:direction=(pos[j]-pos[i])*(1 if voltage[j]>voltage[i] else -1)
    key=(min(i,j),max(i,j),level)
    if key not in cache:
     t=(level-voltage[i])/(voltage[j]-voltage[i]);cache[key]=len(vertices);vertices.append(((1-t)*pos[i]+t*pos[j]).tolist())
    polygon.append(cache[key])
  if len(polygon)<3:continue
  points=np.array([vertices[q] for q in polygon]);center=points.mean(axis=0);normal=np.cross(points[1]-points[0],points[2]-points[0]);length=np.linalg.norm(normal)
  if length<1e-14:continue
  normal/=length
  if normal@direction<0:normal=-normal
  u=points[0]-center;u/=np.linalg.norm(u);v=np.cross(normal,u);order=np.argsort(np.arctan2((points-center)@v,(points-center)@u));faces.append([polygon[q] for q in order]);value=np.mean(dose[ids]);materials.append(min(15,int(15*np.clip(value/max(dose.max(),1e-9),0,1)**.45)))
me=bpy.data.meshes.new('Computed voltage contours');me.from_pydata(vertices,[],faces);me.update();o=bpy.data.objects.new('Potential in surviving 3D material',me);bpy.context.collection.objects.link(o)
for m in mats:me.materials.append(m)
for f,index in zip(me.polygons,materials):f.material_index=index;f.use_smooth=True
# Finite optical foil thickness is styling; contour locations are computed voltage.
mod=o.modifiers.new('Optical observation thickness','SOLIDIFY');mod.thickness=.004
# Actual exposed tetrahedral boundary: no smooth enclosing box or invented shell.
fm=[v for v in boundary.values() if v is not None];glass,gs=material('Smoked violet observation medium',(.31,.23,.49,1),.22,0,.90,1.34)
me=bpy.data.meshes.new('Actual fractured material boundary');me.from_pydata(pos.tolist(),[],fm);me.update();bm=bmesh.new();bm.from_mesh(me);bmesh.ops.recalc_face_normals(bm,faces=list(bm.faces));bmesh.ops.bisect_plane(bm,geom=list(bm.verts)+list(bm.edges)+list(bm.faces),dist=1e-6,plane_co=(0,.18,0),plane_no=(0,1,0),clear_inner=True,clear_outer=False);bm.to_mesh(me);bm.free();o=bpy.data.objects.new('Open optical back section of fractured material',me);bpy.context.collection.objects.link(o);me.materials.append(glass)
sc=bpy.context.scene;sc.render.engine='CYCLES';sc.render.threads_mode='FIXED';sc.render.threads=2;sc.cycles.samples=a.samples;sc.cycles.use_denoising=True;sc.cycles.max_bounces=12;sc.cycles.transmission_bounces=8
if a.device=='METAL':
 prefs=bpy.context.preferences.addons['cycles'].preferences;prefs.compute_device_type='METAL';prefs.get_devices()
 for dev in prefs.devices:dev.use=dev.type=='METAL'
 sc.cycles.device='GPU'
sc.render.resolution_x=a.size;sc.render.resolution_y=a.size;sc.render.resolution_percentage=100;sc.world.use_nodes=True;sc.world.node_tree.nodes['Background'].inputs[0].default_value=(.10,.12,.15,1);sc.world.node_tree.nodes['Background'].inputs[1].default_value=.3;sc.view_settings.view_transform='AgX';sc.view_settings.look='AgX - Medium High Contrast'
def aim(o):o.rotation_euler=(-o.location).to_track_quat('-Z','Y').to_euler()
bpy.ops.object.camera_add(location=(2.8,-4,2.6));cam=bpy.context.object;cam.data.type='ORTHO';cam.data.ortho_scale=2.85;aim(cam);sc.camera=cam
for name,xyz,energy,color,size in [('Warm textile',(-2,-3,4),800,(1,.83,.66),3),('Violet edge',(3,1,2),900,(.57,.65,1),2),('Silver fill',(-3,2,0),450,(.8,.9,1),2)]:
 bpy.ops.object.light_add(type='AREA',location=xyz);o=bpy.context.object;o.name=name;o.data.energy=energy;o.data.color=color;o.data.size=size;aim(o)
a.output.parent.mkdir(parents=True,exist_ok=True);sc.render.image_settings.file_format='PNG';sc.render.filepath=str(a.output.resolve());bpy.ops.render.render(write_still=True)
a.output.with_suffix('.optical.json').write_text(json.dumps({'size':a.size,'samples':a.samples,'device':a.device,'source_nodes':len(pos),'source_tetrahedra_surviving':len(tets),'source_tetrahedra_removed':removed,'contour_faces':len(faces),'levels':isos,'interpretation':'Computed voltage iso-surfaces inside surviving six-edge tetrahedra, actual deformed coordinates. Foil thickness.004 and dose color are optical staging, not added material. Violet back observation surface is the actual surviving tetrahedral boundary exactly plane-clipped at y=.18.','renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
