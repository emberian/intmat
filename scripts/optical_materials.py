"""Declared display materials with distinct optical behavior, not measured simulation properties.

Pure-data presets can be inspected without Blender. make_material imports bpy lazily.
Procedural patina/film varies shading only; it never displaces simulation geometry.
"""
import copy
PRESETS = {
 'amber-resin': {'label':'Amber resin','inputs':{'Base Color':(1,.9,.68,1),'Roughness':.065,'IOR':1.54,'Transmission Weight':.96,'Coat Weight':.3,'Coat Roughness':.04},'absorption':{'color':(.88,.34,.035,1),'density':.45}},
 'smoky-quartz': {'label':'Smoky quartz','inputs':{'Base Color':(.67,.64,.59,1),'Roughness':.085,'IOR':1.544,'Transmission Weight':.98,'Coat Weight':.12,'Coat Roughness':.025},'absorption':{'color':(.44,.42,.37,1),'density':.7}},
 'beeswax': {'label':'Warm translucent wax','inputs':{'Base Color':(.88,.62,.22,1),'Roughness':.42,'IOR':1.44,'Transmission Weight':.025,'Subsurface Weight':.42,'Subsurface Radius':(1,.5,.18),'Subsurface Scale':.07,'Coat Weight':.08,'Coat Roughness':.3}},
 'frosted-ice': {'label':'Frosted ice','inputs':{'Base Color':(.7,.87,.96,1),'Roughness':.3,'IOR':1.31,'Transmission Weight':.72,'Subsurface Weight':.05,'Subsurface Scale':.025},'absorption':{'color':(.66,.85,.94,1),'density':.12}},
 'iridescent-glass': {'label':'Thin-film iridescent glass','inputs':{'Base Color':(.93,.97,1,1),'Roughness':.035,'IOR':1.46,'Transmission Weight':1.0,'Coat Weight':.22,'Coat Roughness':.025,'Thin Film Thickness':420.0,'Thin Film IOR':1.38},'film':{'noise_scale':2.3,'minimum_nm':240.0,'maximum_nm':680.0}},
 'oxidized-copper': {'label':'Oxidized copper','inputs':{'Base Color':(.1,.36,.25,1),'Metallic':.88,'Roughness':.46,'IOR':1.5,'Transmission Weight':0,'Coat Weight':.06},'patina':{'noise_scale':3.7,'metal':(.55,.22,.065,1),'oxide':(.035,.27,.18,1),'minimum_roughness':.24,'maximum_roughness':.65}},
 'brushed-bronze': {'label':'Brushed bronze','inputs':{'Base Color':(.52,.28,.085,1),'Metallic':.95,'Roughness':.28,'Anisotropic':.55,'IOR':1.5,'Transmission Weight':0,'Coat Weight':.12,'Coat Roughness':.12}},
 'porcelain': {'label':'Glazed porcelain','inputs':{'Base Color':(.9,.85,.73,1),'Roughness':.2,'IOR':1.5,'Transmission Weight':0,'Subsurface Weight':.035,'Subsurface Scale':.018,'Coat Weight':.55,'Coat Roughness':.1}},
 'basalt': {'label':'Matte basalt','inputs':{'Base Color':(.045,.052,.06,1),'Metallic':.035,'Roughness':.76,'IOR':1.52,'Transmission Weight':0,'Coat Weight':0}},
}
def specification(name):
 return copy.deepcopy(PRESETS[name])
def make_material(name,preset):
 import bpy
 spec=specification(preset); material=bpy.data.materials.new(name+' / '+spec['label']);material.use_nodes=True
 nodes=material.node_tree.nodes; links=material.node_tree.links; shader=nodes.get('Principled BSDF'); output=nodes.get('Material Output')
 for key,value in spec['inputs'].items():
  if key not in shader.inputs:raise ValueError('Material input unsupported by this Blender: '+key)
  shader.inputs[key].default_value=value
 coordinates=None
 def noise(scale):
  nonlocal coordinates
  if coordinates is None:coordinates=nodes.new('ShaderNodeTexCoord')
  n=nodes.new('ShaderNodeTexNoise');n.inputs['Scale'].default_value=scale;n.inputs['Detail'].default_value=3
  links.new(coordinates.outputs['Generated'],n.inputs['Vector']);return n.outputs['Fac']
 if 'absorption' in spec:
  a=nodes.new('ShaderNodeVolumeAbsorption');a.inputs['Color'].default_value=spec['absorption']['color'];a.inputs['Density'].default_value=spec['absorption']['density'];links.new(a.outputs['Volume'],output.inputs['Volume'])
 if 'film' in spec:
  f=spec['film'];mapping=nodes.new('ShaderNodeMapRange');mapping.inputs['To Min'].default_value=f['minimum_nm'];mapping.inputs['To Max'].default_value=f['maximum_nm'];links.new(noise(f['noise_scale']),mapping.inputs['Value']);links.new(mapping.outputs['Result'],shader.inputs['Thin Film Thickness'])
 if 'patina' in spec:
  p=spec['patina'];fac=noise(p['noise_scale']);ramp=nodes.new('ShaderNodeValToRGB');ramp.color_ramp.elements[0].position=.32;ramp.color_ramp.elements[0].color=p['metal'];ramp.color_ramp.elements[1].position=.6;ramp.color_ramp.elements[1].color=p['oxide'];links.new(fac,ramp.inputs['Fac']);links.new(ramp.outputs['Color'],shader.inputs['Base Color']);rough=nodes.new('ShaderNodeMapRange');rough.inputs['To Min'].default_value=p['minimum_roughness'];rough.inputs['To Max'].default_value=p['maximum_roughness'];links.new(fac,rough.inputs['Value']);links.new(rough.outputs['Result'],shader.inputs['Roughness'])
 return material,spec
