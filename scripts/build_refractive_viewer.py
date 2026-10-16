#!/usr/bin/env python3
"""Package actual mineral captures and Julia ray paths in the shared volume viewer.
Requires numpy for block-mean display reduction; never advances a simulation.
"""
import argparse,base64,gzip,hashlib,json,pathlib,re,subprocess,sys
import numpy as np
p=argparse.ArgumentParser();p.add_argument('run',type=pathlib.Path);p.add_argument('--tag',default='guided-final');a=p.parse_args();root=a.run
run=json.loads((root/'run.json').read_text());n=run['shape'][0];assert n%64==0
stage=root/'viewer-package';stage.mkdir(exist_ok=True);frames=[];observations=[];allpaths=[]
for frame in run['frames']:
    raw=root/frame['directory']/'solid.f32';solid=np.fromfile(raw,dtype='<f4').reshape((n,n,n),order='F');factor=n//64
    reduced=solid.reshape(64,factor,64,factor,64,factor).mean(axis=(1,3,5),dtype=np.float64).astype('<f4')
    name=f"frame-{frame['step']:06d}";(stage/(name+'-deposit.f32')).write_bytes(reduced.tobytes(order='F'))
    frames.append(dict(name=name,time=frame['time'],solid=float(np.sum(solid,dtype=np.float64))))
    path=root/'optics'/f"{frame['step']:06d}-{a.tag}"/'observation.json';obs=json.loads(path.read_text());assert hashlib.sha256(raw.read_bytes()).hexdigest()==obs['source_solid_sha256']
    rays=[r for r in obs['rays'] if r['wavelength']==550]
    selected=[]
    for ray in rays:
        points=np.asarray(ray['path'],dtype=np.float64);points=np.concatenate((points[::8],points[-1:])) if (len(points)-1)%8 else points[::8]
        selected.append(points);allpaths.append(points)
    observations.append((obs,selected,rays,hashlib.sha256(path.read_bytes()).hexdigest()))
maxabs=max(float(np.max(np.abs(p))) for p in allpaths);quantum=maxabs/32760
rayframes=[]
for obs,paths,rays,digest in observations:
    coords=np.concatenate(paths);packed=np.rint(coords/quantum).astype('<i2');assert np.max(np.abs(packed.astype(float)*quantum-coords))<=quantum*.500001
    screen=next(s for s in obs['screens'] if s['wavelength']==550);counts=np.asarray(screen['counts']);assert counts.max()<65536
    check=next(c for c in obs['checks'] if c['wavelength']==550)
    assert int(counts.sum())+check['outside']+check['unresolved']==check['total']
    reference=np.asarray([[r['launch'],r['constant_index_position']] for r in rays],dtype=np.float32).reshape(-1,3)
    rayframes.append(dict(positions=base64.b64encode(packed.tobytes()).decode(),lengths=[len(p) for p in paths],reference=base64.b64encode(reference.astype('<f4').tobytes()).decode(),screen=base64.b64encode(counts.astype('<u2').tobytes()).decode(),screen_size=obs['screen_pixels'],screen_span=obs['screen_span'],check=check,observation_sha256=digest))
meta=dict(title=run['title'],shape=[64]*3,backend='Metal',native_shape=run['shape'],frames=frames,config=dict(simulation=dict(seed=618,dt=run['dt'])),model='Five recorded mineral states, one fixed set of launch apertures. Compare bent paths with their constant-index reference; turn the material or open its interior.',source_sha256=run['source_sha256'],source_run_sha256=hashlib.sha256((root/'run.json').read_bytes()).hexdigest())
(stage/'run.json').write_text(json.dumps(meta));subprocess.run([sys.executable,str(pathlib.Path(__file__).with_name('build_viewer.py')),str(stage),'--fields','deposit'],check=True)
s=(stage/'viewer.html').read_text();match=re.search(r'(<script id="recording" type="application/json">)(.*?)(</script>)',s);data=json.loads(match[2])
data['fields']=json.loads(gzip.decompress(base64.b64decode(data.pop('fields_gzip')))) if 'fields_gzip' in data else data['fields']
for channel in ('activity','nutrient'):data['fields'][channel]=[None]*len(frames)
data['run']['source_channels']={'deposit':'solid','activity':None,'nutrient':None};data['ray_frames_gzip']=base64.b64encode(gzip.compress(json.dumps(rayframes,separators=(',',':')).encode(),mtime=0)).decode();data['ray_quantum']=quantum;data['ray_encoding']='Actual550nm paths: every eighth saved vertex plus endpoints; signed16-bit positions with shared quantum; maximum per-axis rounding error quantum/2. Constant-index endpoints are recorded Julia controls. Receiver bins are losslessuint16.'
data['fields_gzip']=base64.b64encode(gzip.compress(json.dumps(data.pop('fields'),separators=(',',':')).encode(),mtime=0)).decode()
s=s[:match.start(2)]+json.dumps(data,separators=(',',':'))+s[match.end(2):]
def replace(old,new):
    global s
    assert old in s,old[:100]
    s=s.replace(old,new)
replace("if(data.fields_gzip){","if(data.ray_frames_gzip){const stream=new Blob([Uint8Array.from(atob(data.ray_frames_gzip),c=>c.charCodeAt(0))]).stream().pipeThrough(new DecompressionStream('gzip'));data.ray_frames=JSON.parse(await new Response(stream).text());delete data.ray_frames_gzip;}\nif(data.fields_gzip){")
replace("let b=atob(data.fields[c][n]);decoded[key]=Uint8Array.from(b,x=>x.charCodeAt(0));","let encoded=data.fields[c][n];decoded[key]=encoded===null?new Uint8Array(data.run.shape.reduce((a,b)=>a*b,1)):Uint8Array.from(atob(encoded),x=>x.charCodeAt(0));")
replace('RECORDED JULIA · METAL · 64³',f'NATIVE {n}³ · 5 ACTUAL STATES · JULIA RAYS')
replace('<div class="eyebrow">Active materials / specimen 618</div>','<a href="../../#collection-20261016/refractive-memory">← Collection</a><div class="eyebrow">A mineral changes the path of light</div>')
replace('AMBER MEMORY · CYAN EXCITATION','AMBER MINERAL · CYAN TRACED PATHS')
replace('Volume / amber + cyan','Transparent mineral');replace('Slice / deposit','Slice / solid')
replace('<option value="2">Slice / activity</option><option value="3">Slice / nutrient</option>','')
replace('Slice depth','Cut / slice depth');replace('max="0.95" step="0.025" value="0"','max="1" step="0.025" value="1"')
replace('if(dot(p/vec3(.86,1.,.78),p/vec3(.86,1.,.78))>1.)continue;','if(p.z>slice)continue;')
replace('The ellipsoid is a viewing crop. The crop introduces no physical wall; domain boundaries belong to the recorded model. Slice colors use fixed scales across all frames.','Three apertures stay fixed while the mineral changes. Cyan paths were traced through each frozen index field in Julia; pale reference paths use index1. Turning, clipping and exposure change the view only. The receiving plane counts equal-weight550nm rays; rays outside its window remain in the accounting.')
replace('Recorded field dynamics. Playback uses real Julia snapshots, globally quantized to 8 bits per field. Rendering controls affect light only.','The photograph uses native mineral contours. This transparent context averages the native field to64³, with one fixed8-bit scale. Saved paths are thinned and quantized for display; no motion is invented between the five recorded times. Browser light differs from the photograph’s Cycles optics.')
replace('<strong>0</strong>particles','<strong>550 nm</strong>traced wavelength')
replace('Deposit / persistent material','Solid / retained mineral');replace('Activity / excitable fronts','Paths / fixed launch apertures')
replace("let angles=[.43,-.18]","let angles=[.55,1.12]")
replace('index=4,dirty=true','index=2,dirty=true')
replace('t-last>550','t-last>1200')
replace('<div class="below">','''<div class="settings" style="margin:20px 0"><label>Paths<select id="raymode"><option value="bent">Bent through mineral</option><option value="both">Bent + straight reference</option><option value="reference">Straight reference</option><option value="none">Material only</option></select></label><label>Displayed bundle<select id="raystride"><option value="2">Selected threads</option><option value="1">Every recorded thread</option><option value="3">Sparse threads</option></select></label><label><input id="cliprays" type="checkbox" checked style="width:auto;display:inline"> Reveal paths inside mineral</label></div><div class="below">''')
replace('</footer></main>','''</footer><style>.receiver-panel{display:grid;grid-template-columns:minmax(180px,330px) 1fr;gap:24px;margin-top:24px;align-items:center}@media(max-width:650px){.receiver-panel{grid-template-columns:1fr}.receiver-panel canvas{max-width:330px}}</style><section class="receiver-panel"><canvas id="receiver" width="256" height="256" style="width:100%;height:auto;min-height:0;border:1px solid #23353e;border-radius:12px;cursor:default"></canvas><div><h2 style="font:30px Georgia;color:#f6ead4">Where the rays arrive</h2><p id="receiver-count"></p><p>Fixed receiver at z=1.05; 2.4×2.4 window. Pixel brightness shows ray concentration with a fixed display transfer. This is geometric ray density, not interference or a wave-optics caustic.</p></div></section></main>''')
replace("$('time').value=n;dirty=true;}","$('time').value=n;loadRays(n);dirty=true;}")
replace('gl.drawArrays(gl.TRIANGLES,0,6);dirty=false;','gl.drawArrays(gl.TRIANGLES,0,6);drawRays();dirty=false;')
overlay=r'''
const rayVS=`#version 300 es
precision highp float;in vec3 point;out vec3 world;uniform vec2 angles,resolution;uniform float zoom;
mat3 rot(){float a=angles.x,b=angles.y;return mat3(cos(a),0,-sin(a),0,1,0,sin(a),0,cos(a))*mat3(1,0,0,0,cos(b),sin(b),0,-sin(b),cos(b));}
void main(){world=point;vec3 c=transpose(rot())*point-vec3(0,0,zoom);gl_Position=vec4(c.x*1.8/.63/(resolution.x/resolution.y),c.y*1.8/.63,0.,-c.z);}`;
const rayFS=`#version 300 es
precision highp float;precision highp sampler3D;in vec3 world;out vec4 color;uniform sampler3D dep;uniform vec2 scale;uniform float cut,inside;uniform vec4 tint;
void main(){if(world.z>cut)discard;if(inside>.5){if(any(greaterThan(abs(world),vec3(1))))discard;float d=mix(scale.x,scale.y,texture(dep,world*.5+.5).r);if(d<.045)discard;}color=tint;}`;
let rp=gl.createProgram();gl.attachShader(rp,shader(gl.VERTEX_SHADER,rayVS));gl.attachShader(rp,shader(gl.FRAGMENT_SHADER,rayFS));gl.linkProgram(rp);if(!gl.getProgramParameter(rp,gl.LINK_STATUS))throw Error(gl.getProgramInfoLog(rp));
let rayVAO=gl.createVertexArray(),rayBuffer=gl.createBuffer(),refVAO=gl.createVertexArray(),refBuffer=gl.createBuffer();
function initLine(vao,buf){gl.bindVertexArray(vao);gl.bindBuffer(gl.ARRAY_BUFFER,buf);let pos=gl.getAttribLocation(rp,'point');gl.enableVertexAttribArray(pos);gl.vertexAttribPointer(pos,3,gl.FLOAT,false,0,0);gl.bindVertexArray(null);}
initLine(rayVAO,rayBuffer);initLine(refVAO,refBuffer);let rayLengths=[],rayOffsets=[];
function bytes(text){return Uint8Array.from(atob(text),c=>c.charCodeAt(0));}
function loadRays(n){let f=data.ray_frames[n],raw=bytes(f.positions),dv=new DataView(raw.buffer),points=new Float32Array(raw.length/2);for(let i=0;i<points.length;i++)points[i]=dv.getInt16(2*i,true)*data.ray_quantum;gl.bindBuffer(gl.ARRAY_BUFFER,rayBuffer);gl.bufferData(gl.ARRAY_BUFFER,points,gl.STATIC_DRAW);let refs=bytes(f.reference),rdv=new DataView(refs.buffer),reference=new Float32Array(refs.length/4);for(let i=0;i<reference.length;i++)reference[i]=rdv.getFloat32(4*i,true);gl.bindBuffer(gl.ARRAY_BUFFER,refBuffer);gl.bufferData(gl.ARRAY_BUFFER,reference,gl.STATIC_DRAW);rayLengths=f.lengths;let offset=0;rayOffsets=rayLengths.map(length=>{let out=offset;offset+=length;return out});gl.bindBuffer(gl.ARRAY_BUFFER,buffer);
let counts=bytes(f.screen),cdv=new DataView(counts.buffer),ctx=$('receiver').getContext('2d'),img=ctx.createImageData(f.screen_size,f.screen_size);for(let y=0;y<f.screen_size;y++)for(let x=0;x<f.screen_size;x++){let count=cdv.getUint16(2*(x+f.screen_size*y),true),v=1-Math.exp(-.16*count*f.screen_size*f.screen_size/f.check.total),i=4*(x+(f.screen_size-1-y)*f.screen_size);img.data[i]=Math.round(255*Math.pow(v,.65));img.data[i+1]=Math.round(215*Math.pow(v,.7));img.data[i+2]=Math.round(140*Math.pow(v,.8));img.data[i+3]=255;}ctx.putImageData(img,0,0);$('receiver-count').textContent=`t ${data.run.frames[n].time.toFixed(3)} · ${f.check.binned.toLocaleString()} rays in window · ${f.check.outside.toLocaleString()} outside · ${f.check.unresolved} unresolved of ${f.check.total.toLocaleString()} launched.`;}
function drawRays(){let mode=$('raymode').value;if(mode==='none'||+$('mode').value!==0)return;gl.useProgram(rp);let u=name=>gl.getUniformLocation(rp,name);gl.uniform2fv(u('angles'),angles);gl.uniform2f(u('resolution'),canvas.width,canvas.height);gl.uniform1f(u('zoom'),zoom);gl.uniform1f(u('cut'),+$('slice').value);gl.uniform1i(u('dep'),0);gl.uniform2fv(u('scale'),data.scales.deposit);gl.uniform1f(u('inside'),$('cliprays').checked?1:0);gl.enable(gl.BLEND);gl.blendFunc(gl.SRC_ALPHA,gl.ONE);let stride=+$('raystride').value;
if(mode==='bent'||mode==='both'){gl.bindVertexArray(rayVAO);gl.uniform4f(u('tint'),.08,.65,1.,.75);for(let i=0;i<rayLengths.length;i++)if(i%9%stride===0&&Math.floor(i/9)%stride===0)gl.drawArrays(gl.LINE_STRIP,rayOffsets[i],rayLengths[i]);}
if(mode==='reference'||mode==='both'){gl.bindVertexArray(refVAO);gl.uniform4f(u('tint'),.75,.68,.45,.45);for(let i=0;i<rayLengths.length;i++)if(i%9%stride===0&&Math.floor(i/9)%stride===0)gl.drawArrays(gl.LINES,i*2,2);}
gl.disable(gl.BLEND);gl.bindVertexArray(null);gl.useProgram(program);gl.bindBuffer(gl.ARRAY_BUFFER,buffer);}
for(let id of ['raymode','raystride','cliprays'])$(id).onchange=()=>dirty=true;
'''
replace('load(index);function tick(t)',overlay+'\nload(index);function tick(t)')
(root/'viewer.html').write_text(s);print(root/'viewer.html',len(s.encode()),'bytes', 'max path rounding per axis',quantum/2)
