#!/usr/bin/env python3
"""A macro photograph of existing, source-pinned optical feedback interiors.
Only camera, lighting and observation materials change; no field is evolved.
"""
import sys,os,json,hashlib,runpy,math,shutil
from pathlib import Path
import bpy,numpy as np,openvdb as vdb
from mathutils import Vector
args=sys.argv[sys.argv.index('--')+1:];root=Path(args[0]);optics=Path(args[1]);output=Path(args[2]);base=Path(__file__).with_name('render_refractive_memory.py')
output.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(__file__,output.with_suffix('.renderer.py'))
obs=json.loads((optics/'observation.json').read_text());n=obs['shape'][0];dosepath=root/'frames'/f"{obs['step']:06d}"/'dose.f32'
look=tuple(float(v) for v in os.environ.get('SHORE_TARGET','0,.04,-.42').split(','));eye=tuple(float(v) for v in os.environ.get('SHORE_CAMERA','.45,-1.05,-.10').split(','));iso=float(os.environ.get('SHORE_ISO','.25'))
def aim(obj,target):obj.rotation_euler=(Vector(target)-obj.location).to_track_quat('-Z','Y').to_euler()
def material(name,color,roughness,metallic=0,transmission=0,ior=1.4):
    m=bpy.data.materials.new(name);m.use_nodes=True;bs=m.node_tree.nodes.get('Principled BSDF');bs.inputs['Base Color'].default_value=(*color,1);bs.inputs['Roughness'].default_value=roughness;bs.inputs['Metallic'].default_value=metallic;bs.inputs['Transmission Weight'].default_value=transmission;bs.inputs['IOR'].default_value=ior;return m,bs

def landscape(scene):
    bpy.app.handlers.render_pre.remove(landscape)
    camera=scene.camera;camera.location=eye;aim(camera,look);camera.rotation_euler.rotate_axis('Z',-.22);camera.data.lens=48;camera.data.dof.use_dof=True;camera.data.dof.focus_distance=(Vector(eye)-Vector(look)).length;camera.data.dof.aperture_fstop=11
    for obj in list(bpy.data.objects):
        if obj.type=='LIGHT':bpy.data.objects.remove(obj,do_unlink=True)
    for pos,power,color,size,sy in [((-.5,-.3,.6),35,(.5,.9,1),.5,1.8),((.4,.65,-.2),120,(1,.68,.32),.25,1.4),((-.7,.4,-.6),45,(.2,.7,.65),.4,.7),((.15,-.8,-.25),16,(.7,.9,.8),1.2,1.4)]:
        bpy.ops.object.light_add(type='AREA',location=pos);obj=bpy.context.object;obj.data.energy=power;obj.data.color=color;obj.data.shape='RECTANGLE';obj.data.size=size;obj.data.size_y=sy;aim(obj,look)
    scene.world.node_tree.nodes['Background'].inputs['Color'].default_value=(.04,.065,.06,1);scene.world.node_tree.nodes['Background'].inputs['Strength'].default_value=.4
    skin=bpy.data.materials.get('Warm refractive mineral contour');bs=skin.node_tree.nodes.get('Principled BSDF');bs.inputs['Base Color'].default_value=(.36,.65,.52,1);bs.inputs['Transmission Weight'].default_value=.92;bs.inputs['Roughness'].default_value=.32;bs.inputs['IOR'].default_value=1.16;bs.inputs['Coat Weight'].default_value=0
    for node in skin.node_tree.nodes:
        if node.bl_idname=='ShaderNodeVolumeAbsorption':node.inputs['Color'].default_value=(.1,.35,.22,1);node.inputs['Density'].default_value=.3
    dose=np.fromfile(dosepath,dtype='<f4').reshape((n,n,n),order='F');grid=vdb.FloatGrid();grid.copyFromArray(np.ascontiguousarray(dose));points,tris,quads=grid.convertToPolygons(isovalue=iso,adaptivity=.005);points=(points+.5)*2/n-1
    mesh=bpy.data.meshes.new('Measured absorbed-power shoreline');mesh.from_pydata(points.tolist(),[],tris.tolist()+quads.tolist());mesh.update();obj=bpy.data.objects.new('Absorbed-power sheets',mesh);bpy.context.collection.objects.link(obj)
    m,bs=material('Ivory amber absorption contour',(.86,.72,.43),.34,0,.55,1.25);bs.inputs['Subsurface Weight'].default_value=.2;bs.inputs['Subsurface Radius'].default_value=(.05,.025,.012);mesh.materials.append(m)
    for poly in mesh.polygons:poly.use_smooth=True
    scene.view_settings.exposure=-.3
    if '--save-scene' in args:bpy.ops.wm.save_as_mainfile(filepath=str(output.with_suffix('.blend').resolve()))
bpy.app.handlers.render_pre.append(landscape);runpy.run_path(str(base),run_name='__main__')
meta=output.with_suffix('.png.optical.json');data=json.loads(meta.read_text());data['base_renderer_sha256']=data['renderer_sha256'];data['renderer_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();data['macro_observation']=dict(camera=eye,target=look,lens_mm=48,camera_roll=-.22,fstop=11,absorbed_power_isovalue=iso,absorbed_power_sha256=hashlib.sha256(dosepath.read_bytes()).hexdigest(),meaning='Macro photograph of native saved mineral and absorbed-power contours; color, surface refraction, area lights and depth of field are observation choices, not new material dynamics or solved luminescence.');meta.write_text(json.dumps(data,indent=2))
