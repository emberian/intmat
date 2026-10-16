"""Photograph the actual fractured lamina topology and its Joule-dose record."""
import bpy,sys,json,math,argparse
from pathlib import Path
from mathutils import Vector
p=argparse.ArgumentParser();p.add_argument('state');p.add_argument('output');p.add_argument('--size',type=int,default=1000);p.add_argument('--samples',type=int,default=48);p.add_argument('--dose-max',type=float,default=0);p.add_argument('--ortho-scale',type=float,default=2.8);a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
s=json.loads(Path(a.state).read_text());nx,ny,nz=s['dims']
for key,width in [('position',3),('reference',3),('bonds',2)]:
 v=s[key]
 if isinstance(v[0],(int,float)):s[key]=[v[i:i+width] for i in range(0,len(v),width)]
pos=s['position'];edges=s['bonds']
# JSON3 writes matrices as nested columns; normalize from explicit shape.
if len(pos[0])!=3: raise ValueError('position shape')
edge={tuple(sorted((i-1,j-1))):e for e,(i,j) in enumerate(edges)}
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
sc=bpy.context.scene;sc.render.engine='CYCLES';sc.cycles.samples=a.samples;sc.cycles.use_denoising=True;sc.cycles.max_bounces=12;sc.cycles.transmission_bounces=8;sc.render.threads_mode='FIXED';sc.render.threads=2
sc.render.resolution_x=a.size;sc.render.resolution_y=a.size;sc.render.resolution_percentage=100;sc.render.image_settings.file_format='PNG';sc.world.use_nodes=True;sc.world.node_tree.nodes['Background'].inputs[0].default_value=(.09,.12,.17,1);sc.world.node_tree.nodes['Background'].inputs[1].default_value=.30
sc.view_settings.view_transform='AgX';sc.view_settings.look='AgX - Medium High Contrast';sc.view_settings.exposure=-.65
mats=[]
for k in range(12):
 t=k/11;m=bpy.data.materials.new('Joule dose '+str(k));m.use_nodes=True;n=m.node_tree.nodes.get('Principled BSDF');n.inputs['Base Color'].default_value=(.008+.20*t,.028+.065*t,.035+.015*t,1);n.inputs['Metallic'].default_value=.58+.3*t;n.inputs['Roughness'].default_value=.29-.12*t;n.inputs['Coat Weight'].default_value=.45;mats.append(m)
# Each observed lamina is a face complex of the deformed spring network.
# Face failure uses its four actual perimeter bonds; no drawn crack mask.
maxdose=max(a.dose_max,max(s['dose']),1e-15);surviving=0;missing=0
for z in range(nz):
 faces=[];mi=[]
 for y in range(ny-1):
  for x in range(nx-1):
   q=x+nx*y+nx*ny*z;ids=[q,q+1,q+1+nx,q+nx];es=[edge[tuple(sorted((ids[i],ids[(i+1)%4])))] for i in range(4)]
   if any(s['damage'][e]>.5 for e in es):missing+=1;continue
   faces.append(ids);d=sum(s['dose'][i] for i in ids)/4;strain=max(s['strain'][e] for e in es);mi.append(min(11,int(11*min(1,(d/maxdose)**.40*.7+strain*.9))));surviving+=1
 mesh=bpy.data.meshes.new('Surviving lamina');mesh.from_pydata(pos,[],faces);mesh.update();o=bpy.data.objects.new('Lamina %d'%z,mesh);bpy.context.collection.objects.link(o)
 for m in mats:mesh.materials.append(m)
 for f,k in zip(mesh.polygons,mi):f.material_index=k
 mod=o.modifiers.new('Optical finite plate thickness','SOLIDIFY');mod.thickness=.008;mod=o.modifiers.new('Polished actual fracture edge','BEVEL');mod.width=.002;mod.segments=2
# Surviving interlayer struts expose mechanical connections across delaminations.
curve=bpy.data.curves.new('Surviving through-thickness ties','CURVE');curve.dimensions='3D';curve.bevel_depth=.0017;curve.bevel_resolution=1
for e,(i,j) in enumerate(edges):
 i-=1;j-=1
 if s['damage'][e]>.5 or abs(s['reference'][i][2]-s['reference'][j][2])<.04:continue
 # Only axial ties, not all diagonal springs, are opticalized.
 if abs(i-j)!=nx*ny:continue
 line=curve.splines.new('POLY');line.points.add(1)
 for pt,v in zip(line.points,(pos[i],pos[j])):pt.co=(*v,1)
o=bpy.data.objects.new('Exposed conductive ties',curve);bpy.context.collection.objects.link(o);curve.materials.append(mats[-1])
def aim(o,at=(0,0,0)):o.rotation_euler=(Vector(at)-o.location).to_track_quat('-Z','Y').to_euler()
bpy.ops.object.camera_add(location=(2.3,-3.2,2.7));cam=bpy.context.object;cam.data.type='ORTHO';cam.data.ortho_scale=a.ortho_scale;aim(cam);sc.camera=cam
for name,loc,power,col,size in [('Warm longbox',(-2,-3,4),650,(1,.79,.48),3),('Blue edge',(1.8,2,3),1000,(.28,.66,1),2),('White slit',(-3,1,1),850,(.9,1,1),2)]:
 bpy.ops.object.light_add(type='AREA',location=loc);o=bpy.context.object;o.name=name;o.data.energy=power;o.data.color=col;o.data.shape='RECTANGLE';o.data.size=size;o.data.size_y=.5;aim(o)
m=bpy.data.materials.new('Slate');m.use_nodes=True;m.node_tree.nodes.get('Principled BSDF').inputs['Base Color'].default_value=(.014,.022,.030,1);m.node_tree.nodes.get('Principled BSDF').inputs['Roughness'].default_value=.7;bpy.ops.mesh.primitive_plane_add(size=200,location=(0,0,-.50));bpy.context.object.data.materials.append(m)
Path(a.output).parent.mkdir(parents=True,exist_ok=True);sc.render.filepath=a.output;bpy.ops.render.render(write_still=True)
Path(a.output).with_suffix('.json').write_text(json.dumps({'surviving_faces':surviving,'missing_faces':missing,'mapping':'Surviving deformed perimeter bonds define lamina faces. Dose and measured strain color the optical conductor; surviving axial ties are explicit. Plate thickness is optical staging.','size':a.size,'samples':a.samples},indent=2))
