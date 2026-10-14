"""Photograph one copied 2D arena observation; never reads or writes a live arena.
Run with Blender -b --python scripts/render_metabolic_art.py -- SNAPSHOT OUTPUT --size 1200 --samples 64.
"""
import argparse, hashlib, json, math, os, sys
import numpy as np
import bpy
from mathutils import Vector
p=argparse.ArgumentParser()
p.add_argument('snapshot'); p.add_argument('output')
p.add_argument('--size',type=int,default=480); p.add_argument('--samples',type=int,default=16)
p.add_argument('--agent-radius',type=float,default=.004)
a=p.parse_args(sys.argv[sys.argv.index('--')+1:])
blob=open(a.snapshot,'rb').read(); state=json.loads(blob)
nx,ny,nz=state['grid']['shape']
if nz != 1: raise ValueError('This photograph explicitly targets 2D arena observations')
def field(name):
    entry=state['fields'][name]
    return np.array(entry['values'],dtype=np.float32).reshape(entry['shape'],order='F')
phi=field('phi'); interface=field('interface'); film=field('film')
resource=field('substrate_a')+field('substrate_b')
film_scale=.001
f=np.clip(np.arcsinh(film/film_scale)/np.arcsinh(.005/film_scale),0,1)
r=np.clip(resource/.7,0,1)
fraction=np.clip((phi+1)/2,0,1)
amber=np.array([.96,.49,.075]); jade=np.array([.085,.60,.36])
color=fraction[...,None]*amber+(1-fraction[...,None])*jade
# Shared film darkens a copper-colored membrane; resources clear local absorption.
color=color*(1-.45*f[...,None])+np.array([.65,.23,.035])*.45*f[...,None]
density=22+12*interface+35*f+10*(1-r)

def packed_image(name,values):
    image=bpy.data.images.new(name,width=nx,height=ny,float_buffer=True)
    image.colorspace_settings.name='Non-Color'
    if values.ndim==2: values=np.repeat(values[...,None],3,axis=2)
    rgba=np.concatenate([values,np.ones((nx,ny,1),dtype=np.float32)],axis=2)
    image.pixels.foreach_set(np.ascontiguousarray(rgba.transpose(1,0,2),dtype=np.float32).ravel())
    image.pack(); return image
color_image=packed_image('Recorded phase and shared film optical absorption',color)
density_image=packed_image('Recorded interface film and resource optical density',density)
bpy.ops.object.select_all(action='SELECT'); bpy.ops.object.delete(use_global=False)
# A chosen optical thickness; no z motion or 3D biological geometry is simulated.
thickness=.08
bpy.ops.mesh.primitive_cube_add(size=2)
slab=bpy.context.object; slab.name='Recorded living material thin optical section'
slab.scale=(1,1,thickness/2)
bpy.ops.object.transform_apply(location=False,rotation=False,scale=True)
bevel=slab.modifiers.new('Polished optical section edge','BEVEL'); bevel.width=.01; bevel.segments=4
slab.modifiers.new('Face normals','WEIGHTED_NORMAL')
mat=bpy.data.materials.new('Amber jade phase and copper shared film'); mat.use_nodes=True
nodes=mat.node_tree.nodes; nodes.clear(); links=mat.node_tree.links
out=nodes.new('ShaderNodeOutputMaterial')
coords=nodes.new('ShaderNodeTexCoord')
for image in (color_image,density_image):
    tex=nodes.new('ShaderNodeTexImage'); tex.image=image; tex.interpolation='Linear'; tex.extension='EXTEND'
    links.new(coords.outputs['Generated'],tex.inputs['Vector'])
    if image is color_image: color_tex=tex
    else: density_tex=tex
absorb=nodes.new('ShaderNodeVolumeAbsorption')
links.new(color_tex.outputs['Color'],absorb.inputs['Color'])
links.new(density_tex.outputs['Color'],absorb.inputs['Density'])
scatter=nodes.new('ShaderNodeVolumeScatter'); scatter.inputs['Color'].default_value=(.86,.90,.73,1)
scatter.inputs['Density'].default_value=.08
add=nodes.new('ShaderNodeAddShader'); links.new(absorb.outputs[0],add.inputs[0]); links.new(scatter.outputs[0],add.inputs[1])
links.new(add.outputs[0],out.inputs['Volume'])
surface=nodes.new('ShaderNodeBsdfPrincipled')
surface.inputs['Base Color'].default_value=(.98,.985,.96,1)
surface.inputs['Transmission Weight'].default_value=1
surface.inputs['IOR'].default_value=1.34
surface.inputs['Roughness'].default_value=.03
links.new(surface.outputs[0],out.inputs['Surface']); slab.data.materials.append(mat)

lineages=sorted({agent['lineage'] for agent in state['agents']})
colors=[(.97,.81,.37,1),(.67,.28,.09,1),(.91,.47,.22,1),(.71,.82,.60,1),(.80,.64,.34,1),(.68,.33,.19,1),(.93,.87,.60,1),(.46,.65,.49,1)]
bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=1,radius=1)
template=bpy.context.object; base=template.data; bpy.data.objects.remove(template,do_unlink=True)
meshes={}
for rank,lineage in enumerate(lineages):
    mesh=base.copy(); m=bpy.data.materials.new('Observed lineage '+str(lineage)); m.use_nodes=True
    m.node_tree.nodes.clear()
    shader=m.node_tree.nodes.new('ShaderNodeBsdfDiffuse')
    shader.inputs['Color'].default_value=colors[rank%len(colors)]
    shader.inputs['Roughness'].default_value=.2
    output_node=m.node_tree.nodes.new('ShaderNodeOutputMaterial'); m.node_tree.links.new(shader.outputs[0],output_node.inputs['Surface'])
    mesh.materials.append(m); meshes[lineage]=mesh
for agent in state['agents']:
    marker=bpy.data.objects.new('Recorded inhabitant '+str(agent['id']),meshes[agent['lineage']])
    bpy.context.collection.objects.link(marker)
    cell=agent['cell_position']; marker.location=(2*cell[0]/nx-1,2*cell[1]/ny-1,0)
    marker.scale=(a.agent_radius,)*3
# The planar optical backlight reveals absorptive interior rather than emissive fog.
bpy.ops.mesh.primitive_plane_add(size=1.98,location=(0,0,-thickness/2-.002))
back=bpy.context.object; back.name='Diffused transmission light'
m=bpy.data.materials.new('Warm transmission panel'); m.use_nodes=True
m.node_tree.nodes.clear(); emit=m.node_tree.nodes.new('ShaderNodeEmission'); emit.inputs['Color'].default_value=(1,.94,.78,1); emit.inputs['Strength'].default_value=.8
out=m.node_tree.nodes.new('ShaderNodeOutputMaterial'); m.node_tree.links.new(emit.outputs[0],out.inputs[0]); back.data.materials.append(m)
bpy.ops.mesh.primitive_plane_add(size=200,location=(0,0,-thickness/2-.015))
ground=bpy.context.object; m=bpy.data.materials.new('Quiet charcoal stage'); m.use_nodes=True
s=m.node_tree.nodes.get('Principled BSDF'); s.inputs['Base Color'].default_value=(.018,.022,.025,1); s.inputs['Roughness'].default_value=.6; ground.data.materials.append(m)
def aim(obj,target=(0,0,0)):
    obj.rotation_euler=(Vector(target)-obj.location).to_track_quat('-Z','Y').to_euler()
bpy.ops.object.camera_add(location=(1.8,-2.2,3.0)); camera=bpy.context.object; aim(camera); camera.data.lens=48
bpy.context.scene.camera=camera
for pos,power,width,height in [((-1,-2,3),110,2,.12),((1.5,1,3),80,1.5,.15)]:
    bpy.ops.object.light_add(type='AREA',location=pos); light=bpy.context.object
    light.data.energy=power; light.data.shape='RECTANGLE'; light.data.size=width; light.data.size_y=height; aim(light)
scene=bpy.context.scene; scene.render.engine='CYCLES'; scene.cycles.samples=a.samples; scene.cycles.use_denoising=True
scene.cycles.max_bounces=12; scene.cycles.transmission_bounces=8; scene.cycles.volume_bounces=2
scene.render.resolution_x=a.size; scene.render.resolution_y=a.size; scene.render.resolution_percentage=100
scene.world.use_nodes=True; bg=scene.world.node_tree.nodes.get('Background'); bg.inputs['Color'].default_value=(.035,.045,.048,1); bg.inputs['Strength'].default_value=.25
scene.view_settings.view_transform='AgX'
output=os.path.abspath(a.output); os.makedirs(os.path.dirname(output),exist_ok=True)
scene.render.filepath=output; scene.render.image_settings.file_format='PNG'
metadata={'title':'A shared film / amber and jade','source_snapshot':os.path.abspath(a.snapshot),'snapshot_sha256':hashlib.sha256(blob).hexdigest(),'simulation_time':state['time'],'simulation_dimension':2,'grid_shape':state['grid']['shape'],'grid_extent':state['grid']['extent'],'captured_agents':len(state['agents']),'total_population':state['metrics']['population'],'lineages_observed':lineages,'optical_thickness':thickness,'observed_agent_radius':a.agent_radius,'agent_geometry':'Point observation glyphs at actual x/y positions, displayed at z=0; no inferred cell shape or added connectivity','optical_mapping':{'phase':'fraction=(phi+1)/2 clipped to[0,1], jade-to-amber absorption colors','film':'asinh(film/.001)/asinh(.005/.001) clipped to[0,1], copper tint and optical density','resources':'substrate_a+substrate_b normalized by .7, locally clears absorption','density':'22+12*interface+35*film_transfer+10*(1-resource_transfer)','thickness':'Chosen shallow optical extrusion; no 3D dynamics','emission':'No simulated-field emission; diffused physical staging backlight .8'},'surface':{'shader':'Principled transmission','ior':1.34,'roughness':.03},'agent_material':'Diffuse lineage-color observation glyphs','camera':{'position':[1.8,-2.2,3.0],'lens_mm':48},'size':a.size,'samples':a.samples}
with open(os.path.splitext(output)[0]+'.optical.json','w') as stream: json.dump(metadata,stream,indent=2)
bpy.ops.wm.save_as_mainfile(filepath=os.path.splitext(output)[0]+'.blend')
bpy.ops.render.render(write_still=True)
