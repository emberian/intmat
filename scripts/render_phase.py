"""Photograph a geometric 0.5 phase interface from the actual exported scalar grid.
Run: Blender -b --python scripts/render_phase.py -- outputs/phase outputs/phase/phase-glass.png
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
crop=((2,30),(7,25),(2,30)); slices=tuple(slice(*v) for v in crop); slab=f[slices]
grid=openvdb.FloatGrid(background=0); grid.copyFromArray(np.ascontiguousarray(np.pad(slab,1))); vertices,triangles,quads=grid.convertToPolygons(isovalue=.5,adaptivity=0)
vertices=(vertices-1+np.array([v[0] for v in crop])+.5)*2/np.array(shape)-1
bpy.ops.object.select_all(action='SELECT'); bpy.ops.object.delete(use_global=False)
mesh=bpy.data.meshes.new('Computed structure 0.5 isosurface'); mesh.from_pydata(vertices.tolist(),[],triangles.tolist()+quads.tolist()); mesh.update()
obj=bpy.data.objects.new('Demixed phase — optically cut slab',mesh); bpy.context.collection.objects.link(obj)
for face in mesh.polygons: face.use_smooth=True
mat=bpy.data.materials.new('Amber phase optical interpretation'); mat.use_nodes=True; n=mat.node_tree.nodes; links=mat.node_tree.links; shader=n.get('Principled BSDF')
shader.inputs['Base Color'].default_value=(.85,.52,.17,1); shader.inputs['Roughness'].default_value=.08; shader.inputs['IOR'].default_value=1.46; shader.inputs['Transmission Weight'].default_value=.94; shader.inputs['Metallic'].default_value=0
absorb=n.new('ShaderNodeVolumeAbsorption'); absorb.inputs['Color'].default_value=(.92,.43,.07,1); absorb.inputs['Density'].default_value=1.3; links.new(absorb.outputs[0],n.get('Material Output').inputs['Volume']); mesh.materials.append(mat)
# Lift the specimen clear of its studio base so the lowest real channels remain visible.
obj.location.z=.95
bpy.ops.mesh.primitive_plane_add(size=200); floor=bpy.context.object
floor_mat=bpy.data.materials.new('Warm charcoal studio'); floor_mat.use_nodes=True; s=floor_mat.node_tree.nodes.get('Principled BSDF'); s.inputs['Base Color'].default_value=(.022,.025,.029,1); s.inputs['Roughness'].default_value=.28; floor.data.materials.append(floor_mat)
def aim(o,target): o.rotation_euler=(Vector(target)-o.location).to_track_quat('-Z','Y').to_euler()
bpy.ops.object.camera_add(location=(3.5,-5.8,3.4)); camera=bpy.context.object; aim(camera,(0,0,.98)); camera.data.type='ORTHO'; camera.data.ortho_scale=3.65; bpy.context.scene.camera=camera
for name,loc,energy,size,color in [('Warm silk key',(-3,-4,5),1000,3,(1,.88,.71)),('Cool edge',(3,1,4),1500,2,(.70,.85,1)),('Transmitted rear',(-1,3,2),1100,2,(1,.80,.47))]:
    bpy.ops.object.light_add(type='AREA',location=loc); light=bpy.context.object; light.name=name; light.data.energy=energy; light.data.size=size; light.data.color=color; aim(light,(0,0,1))
scene=bpy.context.scene; scene.render.engine='CYCLES'; scene.cycles.samples=a.samples; scene.cycles.use_denoising=True; scene.cycles.max_bounces=16; scene.cycles.transmission_bounces=12; scene.cycles.transparent_max_bounces=12
scene.world.color=(.12,.12,.12); scene.render.resolution_x=a.size; scene.render.resolution_y=a.size; scene.render.resolution_percentage=100; scene.view_settings.view_transform='AgX'; scene.view_settings.exposure=-.45; scene.render.image_settings.file_format='PNG'; scene.render.filepath=output
mapping={'source':os.path.join(directory,'structure.f32'),'source_shape':shape,'source_order':'F','isosurface':.5,'smoothing':{'passes':1,'separable_periodic_weights':[.125,.75,.125]},'crop_indices':crop,'closure':'zero padding closes selected phase at rectangular crop; cut faces are display geometry','mesh_vertices':len(vertices),'mesh_polygons':len(triangles)+len(quads),'optical':{'base_color':[.85,.52,.17],'roughness':.08,'ior':1.46,'transmission':.94,'metallic':0,'volume_absorption_density':1.3,'exposure':-.45,'emission':0},'interpretation':'Geometric demixing interface, not a measured material, fracture, or mechanical boundary','settings':vars(a)}
mapping['source_sha256']=hashlib.sha256(open(mapping['source'],'rb').read()).hexdigest()
mapping['camera']={'position':list(camera.location),'target':[0,0,.98],'projection':'orthographic','scale':3.65}
mapping['lighting']=[{'name':o.name,'position':list(o.location),'power_watts':o.data.energy,'size':o.data.size,'color':list(o.data.color)} for o in scene.objects if o.type=='LIGHT']
json.dump(mapping,open(output+'.optical.json','w'),indent=2); bpy.ops.wm.save_as_mainfile(filepath=os.path.splitext(output)[0]+'.blend'); bpy.ops.render.render(write_still=True)
