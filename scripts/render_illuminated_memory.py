#!/usr/bin/env python3
"""Observation styling for absorbed power inside the light-remodeled mineral.
Reuses the established renderer; the additional geometry is a measured dose
contour, not another simulation or a reconstructed optical path.
"""
import sys,os,json,hashlib,runpy
from pathlib import Path
import bpy,numpy as np,openvdb as vdb
args=sys.argv[sys.argv.index('--')+1:];root=Path(args[0]);optics=Path(args[1]);output=Path(args[2]);base=Path(__file__).with_name('render_refractive_memory.py')
obs=json.loads((optics/'observation.json').read_text());n=obs['shape'][0];frame=root/'frames'/f"{obs['step']:06d}";dosepath=frame/'dose.f32'
cut=float(args[args.index('--cut-y')+1]) if '--cut-y' in args else -2.
threshold=float(os.environ.get('DOSE_CONTOUR','.5'))
def material_observation(scene):
    bpy.app.handlers.render_pre.remove(material_observation)
    for curve in bpy.data.curves:
        if not curve.name.startswith('Recorded '):continue
        pieces=[]
        for spline in curve.splines:
            pts=np.array([tuple(point.co) for point in spline.points]);kept=np.flatnonzero(pts[:,1]>=cut)
            pieces.extend(pts[indices] for indices in np.split(kept,np.flatnonzero(np.diff(kept)>1)+1) if len(indices)>1)
        curve.splines.clear()
        for pts in pieces:
            spline=curve.splines.new('POLY');spline.points.add(len(pts)-1);spline.points.foreach_set('co',pts.ravel())
    for light in bpy.data.lights:light.energy*=float(os.environ.get('SCULPTURE_LIGHT_GAIN','.4'))
    for obj in bpy.data.objects:
        if obj.type=='LIGHT' and obj.location.x<0 and obj.location.y<0:obj.data.energy*=float(os.environ.get('KEY_LIGHT_GAIN','.45'))
    for name in ('Warm refractive mineral contour','Dense mineral / plum gold inclusion'):
        m=bpy.data.materials.get(name)
        if m:
            bs=m.node_tree.nodes.get('Principled BSDF')
            if name.startswith('Warm'):
                bs.inputs['Base Color'].default_value=(.55,.65,.94,1);bs.inputs['Roughness'].default_value=.16;bs.inputs['IOR'].default_value=1.18;bs.inputs['Coat Weight'].default_value=0;bs.inputs['Specular IOR Level'].default_value=.2
                transparent=m.node_tree.nodes.new('ShaderNodeBsdfTransparent');mix=m.node_tree.nodes.new('ShaderNodeMixShader');mix.inputs[0].default_value=float(os.environ.get('SKIN_TRANSPARENCY','.68'));m.node_tree.links.new(bs.outputs[0],mix.inputs[1]);m.node_tree.links.new(transparent.outputs[0],mix.inputs[2]);m.node_tree.links.new(mix.outputs[0],m.node_tree.nodes.get('Material Output').inputs['Surface'])
                for node in m.node_tree.nodes:
                    if node.bl_idname=='ShaderNodeVolumeAbsorption':node.inputs['Color'].default_value=(.28,.42,.75,1)
            else:
                bs.inputs['Base Color'].default_value=(.035,.06,.1,1);bs.inputs['Metallic'].default_value=.15;bs.inputs['Roughness'].default_value=.3
    dose=np.fromfile(dosepath,dtype='<f4').reshape((n,n,n),order='F').copy();ys=(np.arange(n)+.5)*2/n-1;dose[:,ys<cut,:]=0
    grid=vdb.FloatGrid();grid.copyFromArray(np.ascontiguousarray(dose));points,tris,quads=grid.convertToPolygons(isovalue=threshold,adaptivity=.02);points=(points+.5)*2/n-1
    mesh=bpy.data.meshes.new('Recorded absorbed-power density');mesh.from_pydata(points.tolist(),[],tris.tolist()+quads.tolist());mesh.update();obj=bpy.data.objects.new('Absorption contour',mesh);bpy.context.collection.objects.link(obj)
    m=bpy.data.materials.new('Absorbed-power contour / observation glyph');m.use_nodes=True;bs=m.node_tree.nodes.get('Principled BSDF');bs.inputs['Base Color'].default_value=(.12,.65,.53,1);bs.inputs['Metallic'].default_value=.25;bs.inputs['Roughness'].default_value=.27;bs.inputs['Emission Color'].default_value=(.06,.4,.28,1);bs.inputs['Emission Strength'].default_value=.18;mesh.materials.append(m)
    for poly in mesh.polygons:poly.use_smooth=True
    if '--save-scene' in args:bpy.ops.wm.save_as_mainfile(filepath=str(output.with_suffix('.blend').resolve()))
bpy.app.handlers.render_pre.append(material_observation)
runpy.run_path(str(base),run_name='__main__')
meta=output.with_suffix('.png.optical.json');d=json.loads(meta.read_text());d['base_renderer_sha256']=d['renderer_sha256'];d['renderer_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();d['absorbed_power_display']=dict(source_sha256=hashlib.sha256(dosepath.read_bytes()).hexdigest(),isovalue=threshold,skin_transparency=float(os.environ.get('SKIN_TRANSPARENCY','.68')),light_gain=float(os.environ.get('SCULPTURE_LIGHT_GAIN','.4')),key_light_gain=float(os.environ.get('KEY_LIGHT_GAIN','.45')),units='arbitrary energy per material volume per time',meaning='Frozen illumination absorbed-power density used by the preceding material interval. Contour color/emission are observation glyphs, not solved luminescence.');meta.write_text(json.dumps(d,indent=2))
