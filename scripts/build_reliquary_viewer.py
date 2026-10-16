#!/usr/bin/env python3
"""Adapt native saved shells to the existing offline viewer, without simulation.
128³ -> 64³ by exact nonoverlapping 2³ arithmetic means; fixed 8-bit scales.
Three actual states remain discrete; no interpolation between times.
"""
import argparse,array,hashlib,json,pathlib,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('run',type=pathlib.Path);a=p.parse_args()
r=json.loads((a.run/'run.json').read_text());assert r['shape']==[128]*3
out=a.run/'viewer-package';out.mkdir(exist_ok=True)
selected=[f for f in r['frames'] if f['step'] in (0,480,960)];assert len(selected)==3
frames=[];hashes={}
for f in selected:
    name=f"frame-{f['step']:05d}";frames.append(dict(name=name,time=f['time'],solid=r['history'][r['frames'].index(f)]['solid_mean']*128**3))
    for source,alias in [('solid','deposit'),('activity','activity'),('pressure','nutrient')]:
        path=a.run/f['directory']/(source+'.f32');raw=path.read_bytes();assert len(raw)==128**3*4
        hashes[f"{f['step']}/{source}"]=hashlib.sha256(raw).hexdigest();v=array.array('f');v.frombytes(raw)
        if sys.byteorder!='little':v.byteswap()
        reduced=array.array('f')
        for z in range(64):
            for y in range(64):
                for x in range(64):
                    i=2*x+128*(2*y+128*2*z)
                    reduced.append(sum(v[i+dx+128*dy+16384*dz] for dz in (0,1) for dy in (0,1) for dx in (0,1))/8)
        if sys.byteorder!='little':reduced.byteswap()
        (out/(name+'-'+alias+'.f32')).write_bytes(reduced.tobytes())
meta=dict(title=r['title'],shape=[64]*3,backend='Metal',frames=frames,
    model='Turn the three nested mineral shells, cut through them, or inspect slices of solid, recurrent activity and pressure. Choose among three actual recorded states.',
    config=dict(simulation=dict(seed=813,dt=.025)),native_shape=[128]*3,
    source_channels=dict(deposit='solid',activity='activity',nutrient='pressure'),source_sha256=r['source_sha256'],capture_sha256=hashes,
    reduction='Nonoverlapping 2×2×2 arithmetic mean from native128³ to64³, then global linear8-bit quantization per field. Discrete actual times; no temporal interpolation.')
(out/'run.json').write_text(json.dumps(meta))
subprocess.run([sys.executable,str(pathlib.Path(__file__).with_name('build_viewer.py')),str(out)],check=True)
s=(out/'viewer.html').read_text()
def replace(old,new):
    global s
    assert old in s,old[:80]
    s=s.replace(old,new)
replace('<div class="eyebrow">Active materials / specimen 813</div>','<a href="../../#collection-20261015/resonant-reliquary">← Collection</a><div class="eyebrow">What the bell kept / inside the recording</div>')
replace('RECORDED JULIA · METAL · 64³','NATIVE 128³ · VIEW 64³ · 3 ACTUAL STATES')
replace('AMBER MEMORY · CYAN EXCITATION','SOLID · ACTIVITY · PRESSURE')
replace('Slice / nutrient','Slice / pressure')
replace('Slice depth','Cut / slice depth')
replace('max="0.95" step="0.025" value="0"','max="1" step="0.025" value="1"')
replace('The ellipsoid is a viewing crop. The crop introduces no physical wall; domain boundaries belong to the recorded model. Slice colors use fixed scales across all frames.','Drag to rotate; scroll to zoom. The cut slider removes material above the selected z plane. Slice colors use fixed scales across all three times. These are recorded fields, not a live simulation.')
replace('<strong>0</strong>particles','<strong>128³</strong>native voxels')
replace('Recorded field dynamics. Playback uses real Julia snapshots, globally quantized to 8 bits per field. Rendering controls affect light only.','The display averages each 2×2×2 block to 64³, then uses fixed 8-bit field scales. Fine fissures soften at this resolution. Browser light is a scalar observation, separate from the photograph’s Cycles optics.')
replace("if(dot(p/vec3(.86,1.,.78),p/vec3(.86,1.,.78))>1.)continue;",'if(p.z>slice)continue;')
replace("let angles=[.43,-.18]","let angles=[.15,1.35]")
# build_viewer has semantic internal roles; expose the actual channel labels.
replace('"source_channels":{"deposit":"deposit","activity":"activity","nutrient":"nutrient"}','"source_channels":{"deposit":"solid","activity":"activity","nutrient":"pressure"}')
replace('Slice / deposit','Slice / solid');replace('Deposit / persistent material','Solid / retained mineral');replace('Activity / excitable fronts','Activity / recurrent channel')
replace('index=2,dirty=true','index=1,dirty=true')
# A direct contour view makes the retained shell geometry inspectable, while
# volume and signed slices remain available for the other recorded fields.
replace('<option value="0">Volume / amber + cyan</option>','<option value="4" selected>Surface / solid</option><option value="0">Volume / amber + cyan</option>')
replace('if(mode>0){','if(mode>0 && mode<4){')
replace("['deposit','deposit','activity','nutrient'][+$('mode').value]","['deposit','deposit','activity','nutrient','deposit'][+$('mode').value]")
replace('float d=max(f.x,0.),e=max(f.y,0.);',"""float d=max(f.x,0.),e=max(f.y,0.);
if(mode==4){if(d<.12)continue; vec3 tp=p*.5+.5;
vec3 grad=vec3(texture(dep,tp+vec3(voxel.x,0,0)).r-texture(dep,tp-vec3(voxel.x,0,0)).r,texture(dep,tp+vec3(0,voxel.y,0)).r-texture(dep,tp-vec3(0,voxel.y,0)).r,texture(dep,tp+vec3(0,0,voxel.z)).r-texture(dep,tp-vec3(0,0,voxel.z)).r);
vec3 normal=-normalize(grad+vec3(1e-8));if(dot(normal,rd)>0.)normal=-normal;
float light=.22+.68*max(0.,dot(normal,normalize(vec3(-.7,-.8,1.))))+.10*max(0.,dot(normal,normalize(vec3(1.,.4,-.2))));
vec3 pigment=mix(vec3(.16,.42,.43),vec3(.78,.64,.40),clamp((f.y+1.)*.5,0.,1.));radiance=pigment*light;trans=vec3(0);break;}
""")
replace('excitation*2.5','excitation*2.5*min(10.*d,1.)')
replace('Optical density <span','Density (volume) <span')
replace('Slice colors use fixed scales across all three times.','The solid contour is 0.12; its tint follows activity. Slice colors use fixed scales across all three times.')
(a.run/'viewer.html').write_text(s)
print(a.run/'viewer.html',len(s.encode()),'bytes')
