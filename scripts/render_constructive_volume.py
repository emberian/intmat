"""Optical cutaways of actual 3-D cells, bonds and constructed food fields.

The cutaway is an observation choice. Cell radii/overlaps remain physical values;
isosurfaces visualize concentration rather than asserting anatomical membranes.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
import time

import bpy
import numpy as np
import openvdb as vdb
from mathutils import Vector

p = argparse.ArgumentParser()
p.add_argument('run', type=Path)
p.add_argument('output', type=Path)
p.add_argument('--steps', default='all')
p.add_argument('--size', type=int, default=1000)
p.add_argument('--samples', type=int, default=48)
p.add_argument('--device', choices=('CPU', 'METAL'), default='CPU')
p.add_argument('--surface', choices=('retained', 'permeability'), default='retained')
p.add_argument('--surface-optics', choices=('mixed', 'glass'), default='mixed')
p.add_argument('--iso', type=float, default=.45)
p.add_argument('--no-cells', action='store_true')
p.add_argument('--no-cutaway', action='store_true')
p.add_argument('--cut-x', type=float, default=-.08)
p.add_argument('--cut-y', type=float, default=-.08)
p.add_argument('--orbit', type=float, default=.62)
p.add_argument('--elevation', type=float, default=.48)
p.add_argument('--distance', type=float, default=5.6)
p.add_argument('--image-name', default='final.png')
p.add_argument('--save-scene', action='store_true')
a = p.parse_args(sys.argv[sys.argv.index('--') + 1:])
a.output.mkdir(parents=True, exist_ok=True)
run_bytes = (a.run / 'run.json').read_bytes()
run_sha = hashlib.sha256(run_bytes).hexdigest()
(a.output / 'source-run.snapshot.json').write_bytes(run_bytes)
run = json.loads(run_bytes)
selected = None if a.steps == 'all' else set(map(int, a.steps.split(',')))
frames = [f for f in run['frames'] if selected is None or f['step'] in selected]
assert frames, 'No matching captured frame'
renderer_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)
scene = bpy.context.scene
scene.render.engine = 'CYCLES'
scene.cycles.samples = a.samples
scene.cycles.use_denoising = True
scene.cycles.use_adaptive_sampling = True
scene.cycles.adaptive_threshold = .025
scene.cycles.adaptive_min_samples = min(16, a.samples)
scene.cycles.max_bounces = 16
scene.cycles.transmission_bounces = 12
scene.cycles.transparent_max_bounces = 64
if a.device == 'METAL':
    prefs = bpy.context.preferences.addons['cycles'].preferences
    prefs.compute_device_type = 'METAL'
    prefs.get_devices()
    for device in prefs.devices:
        device.use = device.type == 'METAL'
    assert any(d.use for d in prefs.devices)
    scene.cycles.device = 'GPU'
scene.render.resolution_x = scene.render.resolution_y = a.size
scene.render.resolution_percentage = 100
scene.render.image_settings.file_format = 'PNG'
scene.world.use_nodes = True
scene.world.node_tree.nodes['Background'].inputs['Color'].default_value = (.09, .11, .14, 1)
scene.world.node_tree.nodes['Background'].inputs['Strength'].default_value = .3
scene.view_settings.view_transform = 'AgX'
scene.view_settings.look = 'AgX - Medium High Contrast'
scene.view_settings.exposure = -.55


def aim(obj):
    obj.rotation_euler = (-obj.location).to_track_quat('-Z', 'Y').to_euler()


bpy.ops.object.camera_add(location=(a.distance * math.cos(a.elevation) * math.sin(a.orbit),
                                   -a.distance * math.cos(a.elevation) * math.cos(a.orbit),
                                   a.distance * math.sin(a.elevation)))
camera = bpy.context.object
camera.data.lens = 55
scene.camera = camera
cam_position = list(camera.location)
aim(camera)
for name, position, power, color, width in [
    ('Long ivory softbox', (-3, -4, 4), 1100, (1, .92, .75), 3),
    ('Cool interior', (3, 2, 3), 1000, (.48, .75, 1), 2),
    ('Cutaway fill', (3, -4, 1), 700, (.8, .9, 1), 2.5),
    ('Rose edge', (-3, 1, .5), 850, (1, .34, .20), 2),
]:
    bpy.ops.object.light_add(type='AREA', location=position)
    light = bpy.context.object
    light.name = name
    light.data.energy = power
    light.data.color = color
    light.data.shape = 'DISK'
    light.data.size = width
    aim(light)


def material(name, color, metal=0, transmission=0, rough=.2, ior=1.4):
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    shader = mat.node_tree.nodes.get('Principled BSDF')
    for key, value in [('Base Color', color), ('Roughness', rough), ('Metallic', metal),
                       ('Transmission Weight', transmission), ('IOR', ior)]:
        shader.inputs[key].default_value = value
    return mat, shader


def cutaway(mat, base):
    """Transparent clipping at the periodic box and optional front-right wedge."""
    nodes, links = mat.node_tree.nodes, mat.node_tree.links
    geometry = nodes.new('ShaderNodeNewGeometry')
    xyz = nodes.new('ShaderNodeSeparateXYZ')
    links.new(geometry.outputs['Position'], xyz.inputs[0])
    conditions = []
    for axis in 'XYZ':
        absolute = nodes.new('ShaderNodeMath'); absolute.operation = 'ABSOLUTE'
        links.new(xyz.outputs[axis], absolute.inputs[0])
        outside = nodes.new('ShaderNodeMath'); outside.operation = 'GREATER_THAN'
        outside.inputs[1].default_value = 1
        links.new(absolute.outputs[0], outside.inputs[0]); conditions.append(outside.outputs[0])
    if not a.no_cutaway:
        right = nodes.new('ShaderNodeMath'); right.operation = 'GREATER_THAN'; right.inputs[1].default_value = a.cut_x
        front = nodes.new('ShaderNodeMath'); front.operation = 'LESS_THAN'; front.inputs[1].default_value = a.cut_y
        links.new(xyz.outputs['X'], right.inputs[0]); links.new(xyz.outputs['Y'], front.inputs[0])
        both = nodes.new('ShaderNodeMath'); both.operation = 'MULTIPLY'
        links.new(right.outputs[0], both.inputs[0]); links.new(front.outputs[0], both.inputs[1]); conditions.append(both.outputs[0])
    condition = conditions[0]
    for other in conditions[1:]:
        maximum = nodes.new('ShaderNodeMath'); maximum.operation = 'MAXIMUM'
        links.new(condition, maximum.inputs[0]); links.new(other, maximum.inputs[1]); condition = maximum.outputs[0]
    transparent = nodes.new('ShaderNodeBsdfTransparent')
    mix = nodes.new('ShaderNodeMixShader')
    links.new(condition, mix.inputs[0]); links.new(base.outputs[0], mix.inputs[1]); links.new(transparent.outputs[0], mix.inputs[2])
    links.new(mix.outputs[0], nodes.get('Material Output').inputs['Surface'])


skin, skin_shader = material('Actual cell boundary, faint glass', (.65, .78, .85, 1), transmission=.8, rough=.13, ior=1.15)
# Edge-weighted transparency reveals the true sphere surfaces in a dense packing.
nodes, links = skin.node_tree.nodes, skin.node_tree.links
facing = nodes.new('ShaderNodeLayerWeight'); facing.inputs['Blend'].default_value = .25
opacity = nodes.new('ShaderNodeMath'); opacity.operation = 'MULTIPLY_ADD'; opacity.inputs[1].default_value = .22; opacity.inputs[2].default_value = .025
links.new(facing.outputs['Facing'], opacity.inputs[0])
clear = nodes.new('ShaderNodeBsdfTransparent'); skin_mix = nodes.new('ShaderNodeMixShader')
links.new(opacity.outputs[0], skin_mix.inputs[0]); links.new(clear.outputs[0], skin_mix.inputs[1]); links.new(skin_shader.outputs[0], skin_mix.inputs[2])
cutaway(skin, skin_mix)
scaffold, scaffold_shader = material('Retained food interface', (.36, .68, .41, 1), metal=.05, transmission=.9, rough=.12, ior=1.16)
attribute = scaffold.node_tree.nodes.new('ShaderNodeAttribute'); attribute.attribute_name = 'food_a_fraction'
ramp = scaffold.node_tree.nodes.new('ShaderNodeValToRGB')
ramp.color_ramp.elements[0].position = .15; ramp.color_ramp.elements[0].color = (.8, .28, .42, 1)
ramp.color_ramp.elements[1].position = .85; ramp.color_ramp.elements[1].color = (.42, .8, .5, 1)
scaffold.node_tree.links.new(attribute.outputs['Fac'], ramp.inputs[0]); scaffold.node_tree.links.new(ramp.outputs['Color'], scaffold_shader.inputs['Base Color'])
scaffold_shader.inputs['Thin Film Thickness'].default_value = 310
scaffold_shader.inputs['Thin Film IOR'].default_value = 1.32
if a.surface_optics == 'mixed':
    ramp.color_ramp.elements[0].color = (.22, .018, .055, 1)
    ramp.color_ramp.elements[1].color = (.035, .18, .075, 1)
    scaffold_shader.inputs['Metallic'].default_value = .30
    transmission = scaffold.node_tree.nodes.new('ShaderNodeMapRange')
    transmission.inputs['To Min'].default_value = .65
    transmission.inputs['To Max'].default_value = .03
    scaffold.node_tree.links.new(attribute.outputs['Fac'], transmission.inputs['Value'])
    scaffold.node_tree.links.new(transmission.outputs['Result'], scaffold_shader.inputs['Transmission Weight'])
    roughness = scaffold.node_tree.nodes.new('ShaderNodeMapRange')
    roughness.inputs['To Min'].default_value = .12
    roughness.inputs['To Max'].default_value = .28
    scaffold.node_tree.links.new(attribute.outputs['Fac'], roughness.inputs['Value'])
    scaffold.node_tree.links.new(roughness.outputs['Result'], scaffold_shader.inputs['Roughness'])
cutaway(scaffold, scaffold_shader)
bond_mat, bond_shader = material('Actual minimum-image bonds', (.33, .24, .08, 1), metal=.75, rough=.25)
cutaway(bond_mat, bond_shader)
background, bg_shader = material('Matte charcoal cove', (.013, .019, .026, 1), rough=1)
bg_shader.inputs['Specular IOR Level'].default_value = 0
rows = [(-200, -1.20), (3, -1.20)] + [(3 + 4 * math.sin(t), 2.80 - 4 * math.cos(t)) for t in np.linspace(0, math.pi/2, 33)[1:]] + [(7, 200)]
mesh = bpy.data.meshes.new('Cove')
mesh.from_pydata([(x, y, z) for y, z in rows for x in (-200, 200)], [], [(2*q, 2*q+1, 2*q+3, 2*q+2) for q in range(len(rows)-1)])
mesh.update(); obj = bpy.data.objects.new('Cove', mesh); bpy.context.collection.objects.link(obj); mesh.materials.append(background)
for polygon in mesh.polygons: polygon.use_smooth = True
bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=3, radius=1)
sphere = bpy.context.object.data
sphere.name = 'Shared unit cell sphere'
for polygon in sphere.polygons: polygon.use_smooth = True
sphere.materials.append(skin)
bpy.data.objects.remove(bpy.context.object, do_unlink=True)


def add_field(field, ratio):
    # Wrap one voxel so the displayed periodic box does not create an artificial
    # concentration jump at its border. The shader clips at the original box.
    grid = vdb.FloatGrid()
    grid.copyFromArray(np.ascontiguousarray(np.pad(field, 1, mode='wrap')))
    points, triangles, quads = grid.convertToPolygons(isovalue=a.iso, adaptivity=.015)
    indices = np.clip(np.rint(points - 1).astype(int), 0, np.array(field.shape) - 1)
    colors = ratio[tuple(indices.T)] if len(points) else np.empty(0)
    points = (points - .5) * 2 / np.array(field.shape) - 1
    mesh = bpy.data.meshes.new('Saved concentration interface')
    mesh.from_pydata(points.tolist(), [], triangles.tolist() + quads.tolist()); mesh.update()
    obj = bpy.data.objects.new('Saved concentration interface', mesh); bpy.context.collection.objects.link(obj); obj['actual_state'] = True
    mesh.materials.append(scaffold)
    for polygon in mesh.polygons: polygon.use_smooth = True
    attr = mesh.attributes.new('food_a_fraction', 'FLOAT', 'POINT'); attr.data.foreach_set('value', colors.astype(np.float32))
    return len(points)


def add_cells(agents):
    copies = 0
    for cell in agents:
        center = np.asarray(cell['position']) / 6
        radius = cell['radius'] / 6
        shifts = [[0] + ([2] if x-radius < -1 else []) + ([-2] if x+radius > 1 else []) for x in center]
        for shift in itertools.product(*shifts):
            obj = bpy.data.objects.new(f"cell {cell['id']}", sphere)
            bpy.context.collection.objects.link(obj); obj.location = center + np.asarray(shift); obj.scale = (radius,) * 3; obj['actual_state'] = True; copies += 1
    return copies


def add_bonds(agents, edges):
    positions = {c['id']: np.asarray(c['position']) / 6 for c in agents}
    data = bpy.data.curves.new('Periodic physical bond segments', 'CURVE'); data.dimensions = '3D'; data.bevel_depth = .0012; data.bevel_resolution = 1
    for i, j in edges:
        if i not in positions or j not in positions: continue
        start = positions[i]; displacement = positions[j] - start; displacement -= 2 * np.rint(displacement / 2)
        end = start + displacement
        shifts = [[0] + ([2] if min(x,y) < -1 else []) + ([-2] if max(x,y) > 1 else []) for x,y in zip(start,end)]
        for shift in itertools.product(*shifts):
            offset = np.array(shift); spline = data.splines.new('POLY'); spline.points.add(1)
            spline.points[0].co = (*list(start + offset), 1); spline.points[1].co = (*list(end + offset), 1)
    obj = bpy.data.objects.new('Actual stable-ID bonds', data); bpy.context.collection.objects.link(obj); data.materials.append(bond_mat); obj['actual_state'] = True


records = []
for frame in frames:
    started = time.time()
    for obj in list(bpy.data.objects):
        if obj.get('actual_state'):
            data, kind = obj.data, obj.type; bpy.data.objects.remove(obj, do_unlink=True)
            if data is not sphere and data.users == 0: (bpy.data.curves if kind == 'CURVE' else bpy.data.meshes).remove(data)
    directory = a.run / frame['directory']
    volume = json.loads((directory / 'volume.json').read_text())
    assert volume['extent'] == [12.0,12.0,12.0] and volume['order'] == 'F'
    geometry = json.loads((directory / 'geometry.json').read_text())
    shape = tuple(volume['shape'])
    fields = {}
    for name in ('retained_a', 'retained_b', 'permeability'):
        raw = (directory / (name + '.f32')).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == volume['sha256'][name], f'Captured field hash differs: {name}'
        fields[name] = np.frombuffer(raw, dtype=np.float32).reshape(shape, order='F')
    retained = fields['retained_a'] + fields['retained_b']; ratio = fields['retained_a'] / np.maximum(retained, 1e-20)
    surface_points = add_field(retained if a.surface == 'retained' else fields['permeability'], ratio)
    copies = 0 if a.no_cells else add_cells(geometry['agents'])
    add_bonds(geometry['agents'], geometry['bonds'])
    output = a.output / (a.image_name if len(frames) == 1 else f"frame-{frame['step']:07d}.png")
    scene.render.filepath = str(output.resolve())
    inputs = [directory / name for name in ('geometry.json','volume.json','retained_a.f32','retained_b.f32','permeability.f32')]
    optical = {'renderer_sha256': renderer_sha, 'settings': {**vars(a), 'run': str(a.run), 'output': str(a.output)},
               'source': {'run_sha256': run_sha, 'step': frame['step'], 'time': frame['time'], 'files': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}},
               'camera': {'fixed': True, 'location': cam_position, 'lens_mm': 55},
               'geometry': {'physical_cells': len(geometry['agents']), 'periodic_cell_instances': copies, 'bonds': len(geometry['bonds']), 'surface_vertices': surface_points},
               'interpretation': 'Actual 3D cell centers and physical radii, including overlaps, with periodic copies clipped to the original domain. Minimum-image stable-ID bonds. Chosen concentration isosurface tinted by retained A/(A+B), not an anatomical membrane. Mixed optics assigns more opaque jade to A-rich material and clearer rose to B-rich material; this optical contrast is artistic. Front-right optical cutaway reveals the interior. Colors, faint glass skin, fibre radius, thin film and lighting are artistic; no forces, anatomy or temporal interpolation are added.'}
    output.with_suffix('.png.optical.json').write_text(json.dumps(optical, indent=2))
    if a.save_scene: bpy.ops.wm.save_as_mainfile(filepath=str(output.with_suffix('.blend').resolve()))
    bpy.ops.render.render(write_still=True)
    record = {'step': frame['step'], 'time': frame['time'], 'image': output.name, 'sha256': hashlib.sha256(output.read_bytes()).hexdigest(), 'render_seconds': time.time() - started}
    records.append(record); (a.output/'frames.json').write_text(json.dumps({'source_run': str(a.run), 'fixed_camera': True, 'frames': records}, indent=2)); print('FRAME_DONE', json.dumps(record), flush=True)
