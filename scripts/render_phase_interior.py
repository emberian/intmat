"""Photograph a geometric phase-composition isosurface from the actual exported scalar grid.
Run: Blender -b --python scripts/render_phase_interior.py -- outputs/art/phase-window/source outputs/art/phase-window/interior.png
"""
import argparse, hashlib, json, os, sys
import numpy as np
import bpy, openvdb
from mathutils import Vector
p=argparse.ArgumentParser(); p.add_argument('directory'); p.add_argument('output'); p.add_argument('--samples',type=int,default=128); p.add_argument('--size',type=int,default=1000)
a=p.parse_args(sys.argv[sys.argv.index('--')+1:]); directory=os.path.abspath(a.directory); output=os.path.abspath(a.output)
manifest=json.load(open(os.path.join(directory,'volume.json'))); shape=tuple(manifest['shape'])
f=np.fromfile(os.path.join(directory,'structure.f32'),dtype='<f4').reshape(shape,order='F')
# One separable optical smoothing pass suppresses voxel faceting without adding detail.
for axis in range(3):
    f=.125*np.roll(f,1,axis)+.75*f+.125*np.roll(f,-1,axis)
# A rectangular, central slab makes both the cut and the labyrinth legible.
crop=((4,28),(9,23),(4,28)); slices=tuple(slice(*v) for v in crop); slab=f[slices]
grid=openvdb.FloatGrid(background=0); grid.copyFromArray(np.ascontiguousarray(np.pad(slab,1))); vertices,triangles,quads=grid.convertToPolygons(isovalue=.82,adaptivity=0)
vertices=(vertices-1+np.array([v[0] for v in crop])+.5)*2/np.array(shape)-1
bpy.ops.object.select_all(action='SELECT'); bpy.ops.object.delete(use_global=False)
mesh=bpy.data.meshes.new('Computed structure 0.82 contour'); mesh.from_pydata(vertices.tolist(),[],triangles.tolist()+quads.tolist()); mesh.update()
obj=bpy.data.objects.new('Phase-rich interior — composition contour',mesh); bpy.context.collection.objects.link(obj)
for face in mesh.polygons: face.use_smooth=True
mat=bpy.data.materials.new('Amber phase optical interpretation'); mat.use_nodes=True; n=mat.node_tree.nodes; links=mat.node_tree.links; shader=n.get('Principled BSDF')
shader.inputs['Base Color'].default_value=(.55,.73,.66,1); shader.inputs['Roughness'].default_value=.27; shader.inputs['IOR'].default_value=1.38; shader.inputs['Transmission Weight'].default_value=.70; shader.inputs['Metallic'].default_value=0
absorb=n.new('ShaderNodeVolumeAbsorption'); absorb.inputs['Color'].default_value=(.92,.43,.07,1); absorb.inputs['Density'].default_value=.1; links.new(absorb.outputs[0],n.get('Material Output').inputs['Volume']); mesh.materials.append(mat)
# Lift the specimen clear of its studio base so the lowest real channels remain visible.
obj.location.z=.82
# Selected phase is enclosed in a virtual transparent amber display matrix.
bpy.ops.mesh.primitive_cube_add(size=2,location=(0,0,.82)); carrier=bpy.context.object; carrier.name='Transparent amber display matrix'; carrier.scale=(.80,.50,.80); bpy.ops.object.transform_apply(location=False,rotation=False,scale=True)
bevel=carrier.modifiers.new('Polished display edges','BEVEL'); bevel.width=.065; bevel.segments=6; carrier.modifiers.new('Weighted display normals','WEIGHTED_NORMAL')
cm=bpy.data.materials.new('Amber containing matrix'); cm.use_nodes=True; cn=cm.node_tree.nodes; cs=cn.get('Principled BSDF'); cs.inputs['Base Color'].default_value=(1,.93,.78,1); cs.inputs['Transmission Weight'].default_value=1; cs.inputs['IOR'].default_value=1.46; cs.inputs['Roughness'].default_value=.045
ca=cn.new('ShaderNodeVolumeAbsorption'); ca.inputs['Color'].default_value=(.94,.58,.18,1); ca.inputs['Density'].default_value=1.0; cm.node_tree.links.new(ca.outputs[0],cn.get('Material Output').inputs['Volume']); carrier.data.materials.append(cm)
bpy.ops.mesh.primitive_plane_add(size=200); floor=bpy.context.object
floor_mat=bpy.data.materials.new('Warm charcoal studio'); floor_mat.use_nodes=True; s=floor_mat.node_tree.nodes.get('Principled BSDF'); s.inputs['Base Color'].default_value=(.022,.025,.029,1); s.inputs['Roughness'].default_value=.42; floor.data.materials.append(floor_mat)
def aim(o,target): o.rotation_euler=(Vector(target)-o.location).to_track_quat('-Z','Y').to_euler()
bpy.ops.object.camera_add(location=(1.9,-6,1.9)); camera=bpy.context.object; aim(camera,(0,0,.86)); camera.data.type='ORTHO'; camera.data.ortho_scale=2.9; bpy.context.scene.camera=camera
for name,loc,energy,size,color in [('Warm silk key',(-3,-4,2),250,3,(1,.92,.80)),('Cool edge',(3,1,2),180,3,(.65,.85,1)),('Transmitted rear',(-1,3,3),600,3,(1,.88,.65))]:
    bpy.ops.object.light_add(type='AREA',location=loc); light=bpy.context.object; light.name=name; light.visible_camera=False; light.data.energy=energy; light.data.size=size; light.data.color=color; aim(light,(0,0,1))
scene=bpy.context.scene; scene.render.engine='CYCLES'; scene.cycles.samples=a.samples; scene.cycles.use_denoising=True; scene.cycles.max_bounces=16; scene.cycles.transmission_bounces=12; scene.cycles.transparent_max_bounces=12
scene.world.color=(.12,.12,.12); scene.render.resolution_x=a.size; scene.render.resolution_y=a.size; scene.render.resolution_percentage=100; scene.view_settings.view_transform='AgX'; scene.view_settings.exposure=-.3; scene.render.image_settings.file_format='PNG'; scene.render.filepath=output
mapping={'source':os.path.join(directory,'structure.f32'),'source_shape':shape,'source_order':'F','isosurface':.82,'smoothing':{'passes':1,'separable_periodic_weights':[.125,.75,.125]},'crop_indices':crop,'closure':'zero padding closes selected phase at rectangular crop; cut faces are display geometry','mesh_vertices':len(vertices),'mesh_polygons':len(triangles)+len(quads),'optical':{'base_color':[.55,.73,.66],'roughness':.27,'ior':1.38,'transmission':.70,'metallic':0,'volume_absorption_density':.1,'exposure':-.3,'carrier':{'virtual_display_matrix':True,'half_dimensions':[.8,.5,.8],'ior':1.46,'transmission':1,'roughness':.045,'absorption':1.0,'absorption_color':[.94,.58,.18]},'emission':0},'interpretation':'Phase-rich structure=0.82 display contour; not the 0.5 demixing boundary. Carrier and optical properties are display choices, not measured material or mechanical claims.','settings':vars(a)}
mapping['source_sha256']=hashlib.sha256(open(mapping['source'],'rb').read()).hexdigest()
mapping['camera']={'position':list(camera.location),'target':[0,0,.86],'projection':'orthographic','scale':2.9}
mapping['lighting']=[{'name':o.name,'position':list(o.location),'power_watts':o.data.energy,'size':o.data.size,'color':list(o.data.color)} for o in scene.objects if o.type=='LIGHT']
json.dump(mapping,open(output+'.optical.json','w'),indent=2); bpy.ops.wm.save_as_mainfile(filepath=os.path.splitext(output)[0]+'.blend'); bpy.ops.render.render(write_still=True)
