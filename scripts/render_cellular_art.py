"""Photograph physical cells and real junctions from an immutable 2D snapshot.

Blender -b --python scripts/render_cellular_art.py -- snapshot.json image.png
Depth and optics are staging; positions, radii, chemistry and bonds are recorded.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import bpy
import numpy as np
from mathutils import Vector

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('snapshot', type=Path)
p.add_argument('output', type=Path)
p.add_argument('--size', type=int, default=2000)
p.add_argument('--samples', type=int, default=96)
p.add_argument('--device', choices=['CPU', 'METAL'], default='METAL')
a = p.parse_args(sys.argv[sys.argv.index('--') + 1:])
blob = a.snapshot.read_bytes()
state = json.loads(blob)
assert state['grid']['shape'][2] == 1, 'This renderer observes 2D cellular state.'
assert all('radius' in cell for cell in state['agents'])
lx, ly, _ = state['grid']['extent']
scale = 2 / lx
hx, hy = 1., ly / lx
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)


def material(name, color, roughness=.25, metallic=0, transmission=0):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    s = m.node_tree.nodes.get('Principled BSDF')
    s.inputs['Base Color'].default_value = (*color, 1)
    s.inputs['Roughness'].default_value = roughness
    s.inputs['Metallic'].default_value = metallic
    s.inputs['Transmission Weight'].default_value = transmission
    s.inputs['IOR'].default_value = 1.38
    return m


def box(name, location, dimensions, mat, bevel=.01):
    bpy.ops.mesh.primitive_cube_add(size=1, location=location)
    o = bpy.context.object
    o.name = name
    o.dimensions = dimensions
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    o.data.materials.append(mat)
    if bevel:
        m = o.modifiers.new('Optical edge finish', 'BEVEL')
        m.width, m.segments = bevel, 3
        o.modifiers.new('Normals', 'WEIGHTED_NORMAL')
    return o


# Exact recorded grid, bilinear texture interpolation only; no procedural detail.
nx, ny = state['fields']['precursor_a']['shape']
cell_volume = lx * ly / (nx * ny)
fields = [np.asarray(state['fields'][key]['values'], dtype=np.float32).reshape((nx, ny), order='F') / cell_volume
          for key in ('precursor_a', 'precursor_b')]
va, vb = fields
contrast = max(float(np.quantile(va + vb, .98)), 1e-8)
light = np.clip((va + vb) / contrast, 0, 1)
mix = va / (va + vb + 1e-8)
amber = np.array([.45, .20, .035])
teal = np.array([.025, .24, .24])
rgb = (mix[..., None] * amber + (1 - mix[..., None]) * teal) * (.25 + .75 * light[..., None])
rgba = np.concatenate([rgb, np.ones((nx, ny, 1))], axis=2)
img = bpy.data.images.new('Recorded complementary precursor concentrations', width=nx, height=ny, float_buffer=True)
img.colorspace_settings.name = 'Non-Color'
img.pixels.foreach_set(np.ascontiguousarray(rgba.transpose(1, 0, 2), dtype=np.float32).ravel())
img.pack()
bed = material('Recorded resource bath', (.1, .1, .1), .32)
nodes, links = bed.node_tree.nodes, bed.node_tree.links
tex = nodes.new('ShaderNodeTexImage'); tex.image = img; tex.interpolation = 'Linear'
coords = nodes.new('ShaderNodeTexCoord')
links.new(coords.outputs['Generated'], tex.inputs['Vector'])
links.new(tex.outputs['Color'], nodes.get('Principled BSDF').inputs['Base Color'])
box('Chemical bath: recorded 2D resource texture', (0, 0, -.046), (2, 2 * hy, .035), bed)
box('Clear optical section', (0, 0, -.014), (2, 2 * hy, .035), material('Clear gel viewing medium', (.91, .98, .96), .09, transmission=.95))

# Each cell has its measured planar radius. The flattened lens height is optical.
bpy.ops.mesh.primitive_uv_sphere_add(segments=32, ring_count=16)
template = bpy.context.object
base = template.data
bpy.data.objects.remove(template, do_unlink=True)
for polygon in base.polygons: polygon.use_smooth = True
palette = []
for i in range(32):
    t = i / 31
    color = (1-t) * np.array([.09, .63, .55]) + t * np.array([.94, .55, .14])
    m = material('A/B uptake allocation %02d' % i, color, .2, .12, .18)
    m.node_tree.nodes.get('Principled BSDF').inputs['Coat Weight'].default_value = .28
    mesh = base.copy(); mesh.materials.append(m); palette.append(mesh)
centers = {}
for cell in state['agents']:
    x, y = cell['position'][:2]
    centers[cell['id']] = (x, y)
    radius = cell['radius'] * scale
    allocations = cell['metabolism']['allocation']
    t = allocations[0] / max(allocations[0] + allocations[1], 1e-9)
    mesh = palette[min(31, max(0, round(t * 31)))]
    # Periodic images appear only when a body overlaps the observational edge.
    for ox in (-lx, 0, lx):
        for oy in (-ly, 0, ly):
            px, py = (x + ox) * scale, (y + oy) * scale
            if abs(px) > hx + radius or abs(py) > hy + radius: continue
            o = bpy.data.objects.new('Cell %s periodic offset %s %s' % (cell['id'], ox, oy), mesh)
            bpy.context.collection.objects.link(o)
            o.location = (px, py, radius * .22 + .005)
            o.scale = (radius, radius, radius * .32)


def clip_segment(p0, p1):
    low, high = 0., 1.
    d = p1 - p0
    for axis, bound in ((0, hx), (1, hy)):
        if abs(d[axis]) < 1e-12:
            if abs(p0[axis]) > bound: return None
        else:
            u, v = sorted(((-bound-p0[axis])/d[axis], (bound-p0[axis])/d[axis]))
            low, high = max(low, u), min(high, v)
    return None if low >= high else (p0 + low*d, p0 + high*d)


bond_material = material('Actual junction: pale conductive filament', (.76, .9, .77), .22, .5)
visible_segments = 0
for first, second in state['metrics']['cell_bonds']:
    if first not in centers or second not in centers: continue
    start = np.array(centers[first]); delta = np.array(centers[second]) - start
    delta -= np.rint(delta / np.array([lx, ly])) * np.array([lx, ly])
    for ox in (-lx, 0, lx):
        for oy in (-ly, 0, ly):
            s = (start + [ox, oy]) * scale
            clipped = clip_segment(s, s + delta * scale)
            if clipped is None: continue
            curve = bpy.data.curves.new('Junction %s-%s' % (first, second), 'CURVE')
            curve.dimensions = '3D'; curve.bevel_depth = .002; curve.bevel_resolution = 3
            spline = curve.splines.new('POLY'); spline.points.add(1)
            for point, xy in zip(spline.points, clipped): point.co = (*xy, .018, 1)
            o = bpy.data.objects.new(curve.name, curve); bpy.context.collection.objects.link(o)
            curve.materials.append(bond_material); visible_segments += 1

box('Charcoal plinth', (0, 0, -.085), (2.06, 2*hy+.06, .05), material('Charcoal ceramic', (.012, .022, .025), .34))
bpy.ops.mesh.primitive_plane_add(size=200, location=(0, 0, -.12))
bpy.context.object.data.materials.append(material('Quiet studio ground', (.022, .032, .035), .55))


def aim(o): o.rotation_euler = (Vector((0, 0, 0)) - o.location).to_track_quat('-Z', 'Y').to_euler()


bpy.ops.object.camera_add(location=(.35, -1.25, 3.25))
camera = bpy.context.object; camera.data.type = 'ORTHO'; camera.data.ortho_scale = 2.45; aim(camera)
scene = bpy.context.scene; scene.camera = camera
for position, energy, color, size in [((-1.5, -1, 3), 220, (1, .82, .61), 2), ((1, 1.5, 2), 170, (.57, .83, 1), 1.4)]:
    bpy.ops.object.light_add(type='AREA', location=position)
    light = bpy.context.object; light.data.energy = energy; light.data.color = color; light.data.size = size; aim(light)
scene.render.engine = 'CYCLES'; scene.cycles.samples = a.samples; scene.cycles.use_denoising = True
scene.cycles.adaptive_threshold = .02
if a.device == 'METAL':
    prefs = bpy.context.preferences.addons['cycles'].preferences
    prefs.compute_device_type = 'METAL'; prefs.get_devices()
    for device in prefs.devices: device.use = device.type == 'METAL'
    scene.cycles.device = 'GPU'
scene.world.use_nodes = True
scene.world.node_tree.nodes.get('Background').inputs['Color'].default_value = (.08, .12, .14, 1)
scene.world.node_tree.nodes.get('Background').inputs['Strength'].default_value = .3
scene.view_settings.view_transform = 'AgX'
scene.render.resolution_x = scene.render.resolution_y = a.size
scene.render.resolution_percentage = 100
a.output = a.output.resolve(); a.output.parent.mkdir(parents=True, exist_ok=True)
scene.render.filepath = str(a.output); scene.render.image_settings.file_format = 'PNG'
metadata = dict(snapshot_sha256=hashlib.sha256(blob).hexdigest(), source=str(a.snapshot.resolve()), step=state['step'], time=state['time'],
                cells=len(state['agents']), bonds=len(state['metrics']['cell_bonds']), drawn_bond_segments=visible_segments,
                geometry='Recorded physical xy positions and radii; periodic boundary copies. Lens height .32 radius is optical staging, not simulated 3D anatomy.',
                colors='Teal to amber encodes relative import-A versus import-B enzyme allocation, not lineage.',
                medium='Recorded precursor concentrations; hue A/B ratio, brightness total divided by its 98th percentile. Bilinear texture interpolation.',
                resource_contrast=contrast, bond_radius=.002, size=a.size, samples=a.samples)
a.output.with_suffix('.optical.json').write_text(json.dumps(metadata, indent=2))
bpy.ops.wm.save_as_mainfile(filepath=str(a.output.with_suffix('.blend')))
bpy.ops.render.render(write_still=True)
