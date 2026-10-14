"""Blender background entry point: blender -b --python scripts/render_blender.py -- EXPORT_DIR OUTPUT.png [--samples 64] [--size 640] [--orbit 0.55].

Consumes Julia export_volume's real scalar grids; requires Blender's openvdb
(or pyopenvdb) and numpy. Converts to VDB, saves scene.blend, renders Cycles.
"""
import argparse
import json
import math
import os
import sys
import numpy as np
import bpy
from mathutils import Vector
try:
    import openvdb as vdb
except ImportError:
    try:
        import pyopenvdb as vdb
    except ImportError as exc:
        raise RuntimeError("This Blender Python needs openvdb or pyopenvdb for grid conversion.") from exc

parser = argparse.ArgumentParser()
parser.add_argument("directory")
parser.add_argument("output")
parser.add_argument("--samples", type=int, default=64)
parser.add_argument("--size", type=int, default=640)
parser.add_argument("--orbit", type=float, default=0.55)
parser.add_argument("--elevation", type=float, default=0.4)
parser.add_argument("--emission", type=float, default=1.2)
parser.add_argument("--emission-field", default="activity", help="Actual exported scalar driving emission")
parser.add_argument("--palette", choices=("amber-cyan", "amber-gold", "glass-ice", "acoustic-pearl"), default="amber-cyan")
parser.add_argument("--absorption", type=float, default=1.5)
parser.add_argument("--scattering", type=float, default=.28)
parser.add_argument("--scattering-field", default="density")
parser.add_argument("--polished-carrier", action="store_true", help="Optical cut-and-polished ellipsoid fitted to the viewing mask")
parser.add_argument("--ior", type=float, default=1.48)
parser.add_argument("--inclusion-surface", help="Actual exported scalar to mesh as a frosted inclusion")
parser.add_argument("--inclusion-isovalue", type=float, default=.25)
parser.add_argument("--surface-roughness", type=float, default=.035)
parser.add_argument("--density-scale", type=float, default=1.0)
parser.add_argument("--light-power", type=float, default=1.0)
parser.add_argument("--z-scale", type=float, default=1.0, help="Declared shallow optical extrusion scale")
parser.add_argument("--population", help="Actual JSON point population with positions and lineages")
parser.add_argument("--population-radius", type=float, default=.015)
parser.add_argument("--conductive-graph", action="store_true", help="Render face adjacency of binary structure as thin optical filaments")
parser.add_argument("--graph-radius", type=float, default=.008)
parser.add_argument("--bond-graph", help="Actual OBJ vertices and line segments, e.g. deformed surviving bonds")
parser.add_argument("--carrier-shape", choices=("ellipsoid", "box"), default="ellipsoid")
parser.add_argument("--activity-threshold", type=float, default=0.0)
parser.add_argument("--activity-power", type=float, default=1.0)
parser.add_argument("--interfaces", action="store_true", help="Artistic activity/deposit gradient transfer; retains simulation data")
parser.add_argument("--smooth", type=int, default=0, help="Optical three-tap smoothing passes before interface derivatives")
args = parser.parse_args(sys.argv[sys.argv.index("--") + 1:])
directory = os.path.abspath(args.directory)
with open(os.path.join(directory, "volume.json")) as stream:
    manifest = json.load(stream)
shape = tuple(manifest["shape"])
arrays = {}
for name in manifest["fields"]:
    arrays[name] = np.fromfile(os.path.join(directory, name + ".f32"), dtype="<f4").reshape(shape, order="F")
if "density" not in arrays:
    raise ValueError("Export a density field (carrier + mineral), plus optional activity, structure, mask.")
optical_mapping = {"settings": vars(args).copy(), "gradient_scales": {}}
if args.emission_field in arrays:
    arrays["activity"] = arrays[args.emission_field].copy()
elif args.emission_field != "activity":
    raise ValueError(f"Emission source {args.emission_field!r} is absent from the export")
mask = np.clip(arrays.get("mask", np.ones(shape)), 0, 1)
if args.interfaces:
    def edge_strength(field, name):
        for _ in range(args.smooth):
            for axis in range(3):
                before = np.take(field, np.maximum(np.arange(shape[axis]) - 1, 0), axis=axis)
                after = np.take(field, np.minimum(np.arange(shape[axis]) + 1, shape[axis] - 1), axis=axis)
                field = .25 * before + .5 * field + .25 * after
        gradient = [np.gradient(field, 2 / shape[d], axis=d) if shape[d] > 1 else np.zeros_like(field) for d in range(3)]
        magnitude = np.sqrt(sum(component * component for component in gradient))
        scale = max(float(np.quantile(magnitude[mask > .5], .98)), 1e-6)
        optical_mapping["gradient_scales"][name] = scale
        return np.minimum(magnitude / scale, 1.5) ** 1.4
    if "activity" in arrays:
        arrays["activity"] = edge_strength(arrays["activity"], "activity")
    if "structure" in arrays:
        arrays["density"] = .025 + .7 * edge_strength(arrays["structure"], "structure")
arrays["density"] *= args.density_scale
grids = []
if args.scattering_field not in arrays:
    raise ValueError("Scattering source absent from exported arrays")
for name in dict.fromkeys(("density", "activity", "structure", args.scattering_field)):
    if name in arrays:
        grid = vdb.FloatGrid()
        grid.name = name
        values = np.maximum(arrays[name], 0)
        if name == "activity":
            values = np.maximum((values - args.activity_threshold) / max(1e-6, 1 - args.activity_threshold), 0) ** args.activity_power
        grid.copyFromArray(np.ascontiguousarray(values * mask))
        grids.append(grid)
output_path = os.path.abspath(args.output)
os.makedirs(os.path.dirname(output_path), exist_ok=True)
vdb_path = os.path.splitext(output_path)[0] + ".vdb"
vdb.write(vdb_path, grids=grids)
optical_mapping["volume_file"] = vdb_path
output_path = os.path.abspath(args.output)
os.makedirs(os.path.dirname(output_path), exist_ok=True)
for mapping_path in (os.path.join(directory, "optical-mapping.json"), output_path + ".optical.json"):
    with open(mapping_path, "w") as stream:
        json.dump(optical_mapping, stream, indent=2)

bpy.ops.object.select_all(action="SELECT")
bpy.ops.object.delete(use_global=False)
bpy.ops.object.volume_import(filepath=vdb_path)
obj = bpy.context.object
obj.name = "Simulated specimen"
obj.scale = tuple((args.z_scale if d == 2 else 1) * 2 / n for d, n in enumerate(shape))
obj.location = tuple((args.z_scale if d == 2 else 1) * (-1 + 1 / n) for d, n in enumerate(shape))
palettes = {
    "amber-cyan": ((.82,.38,.055,1), (1,.64,.26,1), (.025,.72,1,1)),
    "amber-gold": ((.82,.38,.055,1), (1,.78,.45,1), (1,.56,.13,1)),
    "glass-ice": ((.72,.88,.96,1), (.72,.86,1,1), (.56,.85,1,1)),
    "acoustic-pearl": ((.74,.56,.26,1), (.95,.89,.73,1), (.67,.86,1,1)),
}
absorption_color, scattering_color, emission_color = palettes[args.palette]
material = bpy.data.materials.new(args.palette + " optical interpretation")
material.use_nodes = True
nodes = material.node_tree.nodes
nodes.clear()
links = material.node_tree.links
output = nodes.new("ShaderNodeOutputMaterial")
def attribute(name):
    node = nodes.new("ShaderNodeAttribute")
    node.attribute_name = name
    return node.outputs["Fac"]
def multiply(value, scale):
    node = nodes.new("ShaderNodeMath")
    node.operation = "MULTIPLY"
    links.new(value, node.inputs[0])
    node.inputs[1].default_value = scale
    return node.outputs[0]
def add_shader(a, b):
    node = nodes.new("ShaderNodeAddShader")
    links.new(a, node.inputs[0]); links.new(b, node.inputs[1])
    return node.outputs[0]
absorb = nodes.new("ShaderNodeVolumeAbsorption")
absorb.inputs["Color"].default_value = absorption_color
links.new(multiply(attribute("density"), args.absorption), absorb.inputs["Density"])
scatter = nodes.new("ShaderNodeVolumeScatter")
scatter.inputs["Color"].default_value = scattering_color
scatter.inputs["Anisotropy"].default_value = 0.25
links.new(multiply(attribute(args.scattering_field), args.scattering), scatter.inputs["Density"])
emission = nodes.new("ShaderNodeEmission")
emission.inputs["Color"].default_value = emission_color
if "activity" in arrays and not (args.conductive_graph or args.bond_graph):
    links.new(multiply(attribute("activity"), args.emission), emission.inputs["Strength"])
else:
    emission.inputs["Strength"].default_value = 0
links.new(add_shader(add_shader(absorb.outputs[0], scatter.outputs[0]), emission.outputs[0]), output.inputs["Volume"])
obj.data.materials.append(material)
if args.inclusion_surface:
    if args.inclusion_surface not in arrays:
        raise ValueError("Requested inclusion scalar absent from exported arrays")
    surface_grid = vdb.FloatGrid()
    surface_grid.copyFromArray(np.ascontiguousarray(arrays[args.inclusion_surface] * mask))
    points, triangles, quads = surface_grid.convertToPolygons(isovalue=args.inclusion_isovalue, adaptivity=.05)
    points = 2 * (points + .5) / np.array(shape) - 1
    points[:, 2] *= args.z_scale
    mesh = bpy.data.meshes.new("Actual deposited inclusion isosurface")
    mesh.from_pydata(points.tolist(), [], triangles.tolist() + quads.tolist())
    inclusion = bpy.data.objects.new("Baked field inclusion", mesh)
    bpy.context.collection.objects.link(inclusion)
    inclusion_material = bpy.data.materials.new("Pale frosted deposited phase")
    inclusion_material.use_nodes = True
    shader = inclusion_material.node_tree.nodes.get("Principled BSDF")
    shader.inputs["Base Color"].default_value = (.85,.93,.98,1)
    shader.inputs["Roughness"].default_value = .18
    shader.inputs["IOR"].default_value = 1.2
    shader.inputs["Transmission Weight"].default_value = .3
    inclusion.data.materials.append(inclusion_material)
    for polygon in mesh.polygons:
        polygon.use_smooth = True
    optical_mapping["inclusion_mesh"] = {"vertices": len(points), "triangles": len(triangles), "quads": len(quads), "isovalue": args.inclusion_isovalue}

if args.conductive_graph or args.bond_graph:
    graph_edges = []
    if args.bond_graph:
        graph_points = []
        with open(args.bond_graph) as stream:
            for line in stream:
                parts = line.split()
                if parts and parts[0] == "v":
                    graph_points.append(Vector(tuple(float(value) for value in parts[1:4])))
                elif parts and parts[0] == "l":
                    ids = [int(value.split("/")[0]) - 1 for value in parts[1:]]
                    for left, right in zip(ids, ids[1:]):
                        graph_edges.append((graph_points[left], graph_points[right]))
    else:
        if "structure" not in arrays:
            raise ValueError("Conductive graph requires the actual binary conductor structure field")
        occupied = arrays["structure"] > .5
        for index in np.argwhere(occupied):
            for axis in range(3):
                neighbor = index.copy()
                neighbor[axis] += 1
                if neighbor[axis] >= shape[axis] or not occupied[tuple(neighbor)]:
                    continue
                p = Vector(tuple(-1 + 2 * (index[d] + .5) / shape[d] for d in range(3)))
                q = Vector(tuple(-1 + 2 * (neighbor[d] + .5) / shape[d] for d in range(3)))
                graph_edges.append((p, q))
    vertices, faces = [], []
    for p, q in graph_edges:
        direction = (q - p).normalized()
        reference = Vector((0, 1, 0) if abs(direction.y) < .9 else (1, 0, 0))
        u = direction.cross(reference).normalized()
        v = direction.cross(u)
        start = len(vertices)
        for center in (p, q):
            for side in range(8):
                angle = 2 * math.pi * side / 8
                point = center + args.graph_radius * (math.cos(angle) * u + math.sin(angle) * v)
                vertices.append(tuple(point))
        for side in range(8):
            following = (side + 1) % 8
            faces.append((start + side, start + following, start + 8 + following, start + 8 + side))
        faces.extend((tuple(start + side for side in reversed(range(8))), tuple(start + 8 + side for side in range(8))))
    mesh = bpy.data.meshes.new("Actual conductive face adjacency")
    mesh.from_pydata(vertices, [], faces)
    filament = bpy.data.objects.new("Conductive lattice optical filaments", mesh)
    bpy.context.collection.objects.link(filament)
    filament_material = bpy.data.materials.new("Gold conductor optical surface")
    filament_material.diffuse_color = emission_color
    filament_material.use_nodes = True
    shader = filament_material.node_tree.nodes.get("Principled BSDF")
    shader.inputs["Base Color"].default_value = emission_color
    shader.inputs["Metallic"].default_value = .25
    shader.inputs["Roughness"].default_value = .28
    shader.inputs["Emission Color"].default_value = emission_color
    shader.inputs["Emission Strength"].default_value = .2 * args.emission
    filament.data.materials.append(filament_material)
    for polygon in mesh.polygons:
        polygon.use_smooth = True

if args.population:
    with open(args.population) as stream:
        population = json.load(stream)
    lineages = sorted(set(population["lineages"]))
    colors = [(1,.48,.08,1), (.06,.68,.78,1), (.65,.2,.35,1), (.25,.58,.24,1), (.47,.25,.75,1), (.82,.75,.43,1)]
    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=2, radius=1)
    template = bpy.context.object
    base_mesh = template.data
    bpy.data.objects.remove(template, do_unlink=True)
    marker_meshes = {}
    for rank, lineage in enumerate(lineages):
        mesh = base_mesh.copy()
        marker_material = bpy.data.materials.new("Recorded lineage " + str(lineage))
        marker_material.diffuse_color = colors[rank % len(colors)]
        marker_material.use_nodes = True
        shader = marker_material.node_tree.nodes.get("Principled BSDF")
        shader.inputs["Base Color"].default_value = colors[rank % len(colors)]
        shader.inputs["Roughness"].default_value = .28
        shader.inputs["Metallic"].default_value = .2
        mesh.materials.append(marker_material)
        marker_meshes[lineage] = mesh
    for identity, position, lineage in zip(population["ids"], population["positions"], population["lineages"]):
        marker = bpy.data.objects.new("Recorded agent " + str(identity), marker_meshes[lineage])
        bpy.context.collection.objects.link(marker)
        marker.location = (2 * position[0] / shape[0] - 1, 2 * position[1] / shape[1] - 1, 0)
        marker.scale = (args.population_radius,) * 3
    optical_mapping["population_observation"] = {"count": len(population["ids"]), "lineages": lineages, "marker_radius": args.population_radius, "source": args.population, "dimension": 2, "optical_extrusion_z_scale": args.z_scale}

if args.polished_carrier:
    selected = np.argwhere(mask > .01)
    if selected.size == 0:
        raise ValueError("Cannot cut a polished carrier from an empty viewing mask")
    axes = tuple(max(abs(-1 + 2 * (selected[:, d].min() + .5) / shape[d]), abs(-1 + 2 * (selected[:, d].max() + .5) / shape[d])) + 1 / shape[d] for d in range(3))
    if args.carrier_shape == "ellipsoid":
        bpy.ops.mesh.primitive_uv_sphere_add(segments=96, ring_count=64, radius=1)
    else:
        bpy.ops.mesh.primitive_cube_add(size=2)
    shell = bpy.context.object
    shell.name = "Optically cut polished specimen"
    axes = (axes[0], axes[1], axes[2] * args.z_scale)
    shell.scale = axes
    if args.carrier_shape == "box":
        bevel = shell.modifiers.new("Polished edge bevel", "BEVEL")
        bevel.width = .025
        bevel.segments = 4
        shell.modifiers.new("Face normals", "WEIGHTED_NORMAL")
    for polygon in shell.data.polygons:
        polygon.use_smooth = args.carrier_shape == "ellipsoid"
    shell_material = bpy.data.materials.new("Polished refractive optical cut")
    shell_material.use_nodes = True
    shader = shell_material.node_tree.nodes.get("Principled BSDF")
    shader.inputs["Base Color"].default_value = (.98, .975, .96, 1)
    shader.inputs["Roughness"].default_value = args.surface_roughness
    shader.inputs["IOR"].default_value = args.ior
    shader.inputs["Transmission Weight"].default_value = 1
    shell.data.materials.append(shell_material)
    bpy.ops.mesh.primitive_plane_add(size=200, location=(0,0,-axes[2]-.04))
    ground = bpy.context.object
    ground.name = "Dark studio ground"
    ground_material = bpy.data.materials.new("Studio charcoal")
    ground_material.use_nodes = True
    shader = ground_material.node_tree.nodes.get("Principled BSDF")
    shader.inputs["Base Color"].default_value = (.012,.016,.023,1)
    shader.inputs["Roughness"].default_value = .3
    shader.inputs["Metallic"].default_value = .15
    ground.data.materials.append(ground_material)
    optical_mapping["polished_cut_axes"] = list(axes)
    for mapping_path in (os.path.join(directory, "optical-mapping.json"), output_path + ".optical.json"):
        with open(mapping_path, "w") as stream:
            json.dump(optical_mapping, stream, indent=2)

def point_at(obj, target):
    obj.rotation_euler = (Vector(target) - obj.location).to_track_quat("-Z", "Y").to_euler()
bpy.ops.object.camera_add(location=(4.65 * math.cos(args.elevation) * math.sin(args.orbit), -4.65 * math.cos(args.elevation) * math.cos(args.orbit), 4.65 * math.sin(args.elevation)))
camera = bpy.context.object
point_at(camera, (0, 0, 0))
camera.data.lens = 52
bpy.context.scene.camera = camera
for location, energy, size, color in [((-3, -2, 4), 650, 4, (1, .74, .4)), ((2, 2, 3), 950, 3, (.45, .72, 1))]:
    bpy.ops.object.light_add(type="AREA", location=location)
    light = bpy.context.object
    light.data.energy = energy * args.light_power
    light.data.shape = "DISK"
    light.data.size = size
    light.data.color = color
    point_at(light, (0, 0, 0))
scene = bpy.context.scene
scene.render.engine = "CYCLES"
scene.cycles.samples = args.samples
scene.cycles.use_denoising = True
scene.cycles.volume_step_rate = 0.25
scene.cycles.volume_max_steps = 2048
scene.cycles.max_bounces = 12
scene.cycles.transmission_bounces = 8
scene.cycles.volume_bounces = 4
scene.render.resolution_x = args.size
scene.render.resolution_y = args.size
scene.render.resolution_percentage = 100
scene.world.use_nodes = True
background = scene.world.node_tree.nodes.get("Background")
background.inputs["Color"].default_value = (0.008, 0.014, 0.025, 1)
background.inputs["Strength"].default_value = 0.3
scene.render.image_settings.file_format = "PNG"
scene.render.filepath = os.path.abspath(args.output)
os.makedirs(os.path.dirname(scene.render.filepath), exist_ok=True)
for mapping_path in (os.path.join(directory, "optical-mapping.json"), output_path + ".optical.json"):
    with open(mapping_path, "w") as stream:
        json.dump(optical_mapping, stream, indent=2)
obj.data.filepath = "//" + os.path.basename(vdb_path)
bpy.ops.wm.save_as_mainfile(filepath=os.path.splitext(output_path)[0] + ".blend")
bpy.ops.render.render(write_still=True)
