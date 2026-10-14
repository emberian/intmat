#!/usr/bin/env python3
"""Package recorded Julia fields into an offline WebGL2 volume viewer.

Each channel is globally, linearly quantized to 8 bits; scales are embedded so
frames retain comparable values. Julia's x-fastest storage maps directly to
WebGL's x-fastest texture layout. No simulation is run in the browser.
"""
import argparse, array, base64, html as html_escape, json, math, pathlib, sys
try:
    import tomllib
except ImportError:
    tomllib = None
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('directory', nargs='?', default='outputs/choir')
p.add_argument('--fields', help='Up to three recorded field names for amber, cyan, and the third diagnostic slice')
a=p.parse_args(); root=pathlib.Path(a.directory)
metadata=json.loads((root/'run.json').read_text())
streamed='budget' in metadata and (root/'frames.jsonl').exists()
run=dict(metadata)
if streamed:
    budget=metadata['budget']; n=int(budget['resolution'])
    profile_path=root/'profile.toml'
    if not tomllib:
        raise SystemExit('Reading big-run profile.toml requires Python 3.11 or later.')
    with profile_path.open('rb') as profile_file: config=tomllib.load(profile_file)
    frames=[]
    for line_number,line in enumerate((root/'frames.jsonl').read_text().splitlines(),1):
        if not line.strip():continue
        item=json.loads(line)
        frames.append(dict(item,name=item['prefix'],time=item['scheduled_time'],solid=item.get('deposit',0)))
    run.update(shape=[n,n,n],backend=budget['backend'],config=config,frames=frames,
        title=config.get('title',root.name.replace('-',' ').title()),
        model='Recorded phase separation fields.' if budget['kind']=='phase' else 'Recorded excitable, diffusive, and deposited material fields.',
        extent=budget.get('extent',config.get('simulation',{}).get('extent')),
        time_basis='scheduled_time = step × Float64(dt)',recording_schema='streamed bigrun')
    run['config'].setdefault('simulation',{})['dt']=budget['dt']
    available=list(config['simulation']['capture_fields'])
    defaults=['phi','interface','nutrient'] if budget['kind']=='phase' else ['deposit','activity','nutrient']
else:
    available=['deposit','activity','nutrient']; defaults=available
shape=run['shape']
if len(shape)!=3 or any(not isinstance(n,int) or n<1 for n in shape):
    raise SystemExit(f'Viewer requires a positive three-dimensional shape, received {shape}')
frames=run['frames']
if not frames:raise SystemExit('No recorded frames to view')
# The stream order is the captured order. Reject accidental clock regressions.
if any(not math.isfinite(f['time']) for f in frames) or any(b['time']<a['time'] for a,b in zip(frames,frames[1:])):
    raise SystemExit('Recorded frame times must be finite and nondecreasing')
requested=a.fields.split(',') if a.fields else defaults
if not 1<=len(requested)<=3:raise SystemExit('--fields requires one to three comma-separated names')
if a.fields and any(c not in available for c in requested):
    raise SystemExit(f'Requested field is not recorded. Available: {", ".join(available)}')
source_channels={role:(requested[i] if i<len(requested) and requested[i] in available else None)
                 for i,role in enumerate(['deposit','activity','nutrient'])}
count=math.prod(shape); fields={};scales={}
def read_field(frame,channel):
    path=root/(frame['name']+'-'+channel+'.f32')
    if path.stat().st_size!=count*4:raise SystemExit(f'{path}: expected {count*4} bytes for shape {shape}')
    v=array.array('f')
    with path.open('rb') as file:v.fromfile(file,count)
    if sys.byteorder!='little':v.byteswap()
    if not all(math.isfinite(x) for x in v):raise SystemExit(f'{path}: nonfinite field value')
    return v
for role,c in source_channels.items():
    if c is None:
        scales[role]=[0,0]; zero=base64.b64encode(bytes(count)).decode();fields[role]=[zero]*len(frames);continue
    # Two passes bound temporary Float32 memory to one field, not the entire run.
    lo=math.inf;hi=-math.inf
    for f in frames:
        v=read_field(f,c);lo=min(lo,min(v));hi=max(hi,max(v))
    scales[role]=[lo,hi]; encoded=[]
    multiplier=255/(hi-lo) if hi>lo else 0
    for f in frames:
        v=read_field(f,c)
        q=bytes(min(255,max(0,round((x-lo)*multiplier))) for x in v)
        encoded.append(base64.b64encode(q).decode())
    fields[role]=encoded
run['source_channels']=source_channels
run['preview_encoding']='global linear 8-bit quantization per field; original Float32 files unchanged'
data=json.dumps(dict(run=run,scales=scales,fields=fields),separators=(',',':')).replace('</','<\\/')
html=r'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><link rel="icon" href="data:,"><title>Choir — a material remembers</title>
<style>
:root{color-scheme:dark;font:15px system-ui,sans-serif;color:#c3c6c7;background:#070d12}*{box-sizing:border-box}body{margin:0}main{max-width:1500px;margin:auto;padding:36px 5vw}header{display:flex;justify-content:space-between;align-items:center;gap:20px}h1{font:400 clamp(35px,5vw,66px) Georgia,serif;color:#f6ead4;margin:8px 0}h1 em{font-weight:400;color:#c1995d}p{line-height:1.6}.eyebrow{letter-spacing:.22em;font-size:11px;color:#93adaf;text-transform:uppercase}.badge{border:1px solid #32434a;padding:10px 14px;border-radius:30px;font-size:11px;white-space:nowrap}.stage{position:relative;margin:24px 0 0;border:1px solid #23353e;border-radius:18px;overflow:hidden;background:#071016}canvas{display:block;width:100%;height:min(66vh,680px);min-height:350px;touch-action:none;cursor:grab}.stage-note{position:absolute;top:20px;left:24px;font-size:11px;letter-spacing:.13em;pointer-events:none}.hint{position:absolute;bottom:20px;right:24px;color:#819499;font-size:12px;pointer-events:none}.controls{padding:22px 26px;border:1px solid #23353e;border-top:0;border-radius:0 0 18px 18px;margin-top:-15px;padding-top:30px;background:#0b141b}.timeline{display:flex;gap:20px;align-items:center}.timeline input{flex:1}button{background:#173039;color:#e1e8e6;border:1px solid #35545f;padding:10px 18px;border-radius:8px;cursor:pointer}input{accent-color:#b7985b}input[type=range]{width:100%}.time{min-width:115px;font:14px ui-monospace,monospace;color:#e5c58f}.settings{display:grid;grid-template-columns:repeat(4,1fr);gap:25px;margin-top:24px}label{font-size:12px;color:#9bb0b5;display:block}label span{float:right;color:#e1c792}label input{display:block;margin-top:12px}select{background:#15242c;color:#d9e1e0;border:1px solid #35545f;border-radius:5px;padding:7px;width:100%;margin-top:7px}.below{display:grid;grid-template-columns:1.4fr 1fr;gap:50px;margin-top:28px;font-size:13px}.legend{display:flex;gap:24px;flex-wrap:wrap}.dot{width:7px;height:7px;border-radius:50%;display:inline-block;margin-right:8px}.stats{display:flex;gap:25px;font:12px ui-monospace,monospace;color:#839da5}.stats strong{color:#c4d1d1;display:block;margin-bottom:7px;font-size:16px}footer{margin-top:25px;color:#607a85;font-size:11px}#error{color:#e7a57c;padding:30px;display:none}@media(max-width:650px){main{padding:22px 16px}header{display:block}.badge{display:inline-block;margin-top:12px}.settings{grid-template-columns:1fr 1fr;gap:15px}.below{display:block}.timeline{gap:10px}.controls{padding:25px 15px}.time{min-width:85px;font-size:12px}.stats{margin-top:25px}}
</style><main><header><div><div class="eyebrow">Active materials / specimen 017</div><h1>Choir <em>/ a material remembers</em></h1></div><div class="badge">RECORDED JULIA · CPU · 40³</div></header><div class="stage"><canvas id="view" aria-label="Interactive volume of recorded deposit and activity. Drag to rotate."></canvas><div class="stage-note">AMBER MEMORY · CYAN EXCITATION</div><div class="hint">Drag to turn · scroll to zoom</div><div id="error"></div></div><div class="controls"><div class="timeline"><button id="play">▶ Play</button><input id="time" aria-label="Recorded frame" type="range" min="0" max="20" value="10" step="1"><output class="time" id="timestamp"></output></div><div class="settings"><label>Optical density <span id="ov"></span><input id="opacity" type="range" min="0.1" max="5" step="0.05" value="1.1"></label><label>Exposure <span id="ev"></span><input id="exposure" type="range" min="0.3" max="4" step="0.05" value="1.6"></label><label>View<select id="mode"><option value="0">Volume / amber + cyan</option><option value="1">Slice / deposit</option><option value="2">Slice / activity</option><option value="3">Slice / nutrient</option></select></label><label>Slice depth <span id="sv"></span><input id="slice" type="range" min="-0.95" max="0.95" step="0.025" value="0"></label></div></div><div class="below"><div><div class="legend"><div><i class="dot" style="background:#d99e48"></i>Deposit / persistent material</div><div><i class="dot" style="background:#42d9df"></i>Activity / excitable fronts</div></div><p>An excitable medium and recurrent particle policies leave a material trace. Cyan emphasizes spatial changes in activity inside that accumulated amber. Scrub the recording to watch the balance change.</p><p style="color:#79929b">The ellipsoid is a viewing crop. The simulated domain is a periodic cube; the crop introduces no physical wall. Slice colors use fixed scales across all frames.</p></div><div><div class="stats"><div><strong id="solid"></strong>deposit sum</div><div><strong id="frame"></strong>captured frame</div><div><strong>600</strong>particles</div></div><p style="color:#79929b">Fixed recurrent policies; no training or evolution. Playback uses real Julia snapshots, globally quantized to 8 bits per field. Rendering controls affect light only.</p><div id="scale" style="font:11px ui-monospace,monospace;color:#8ba3ad"></div></div></div><footer>21 recorded frames · seed 17 · dt 0.02 · x fastest · offline WebGL2</footer></main><script id="recording" type="application/json">__DATA__</script><script>
const data=JSON.parse(document.getElementById('recording').textContent),$=id=>document.getElementById(id),canvas=$('view'),gl=canvas.getContext('webgl2',{alpha:false,antialias:false,preserveDrawingBuffer:true});
if(!gl){$('error').style.display='block';$('error').textContent='This volume requires WebGL2. Open this file in a browser with hardware acceleration enabled.';}else{try{
const vertex=`#version 300 es
in vec2 pos;out vec2 uv;void main(){uv=pos;gl_Position=vec4(pos,0,1);}`;
const fragment=`#version 300 es
precision highp float;precision highp sampler3D;
in vec2 uv;out vec4 color;uniform sampler3D dep,act,nut;uniform vec2 resolution;uniform vec3 voxel;uniform vec2 angles;uniform vec2 ds,as,ns;uniform float opacity,exposure,zoom,slice;uniform int mode;
mat3 rot(){float a=angles.x,b=angles.y;return mat3(cos(a),0,-sin(a),0,1,0,sin(a),0,cos(a))*mat3(1,0,0,0,cos(b),sin(b),0,-sin(b),cos(b));}
vec3 sampleFields(vec3 p){vec3 t=p*.5+.5;return vec3(mix(ds.x,ds.y,texture(dep,t).r),mix(as.x,as.y,texture(act,t).r),mix(ns.x,ns.y,texture(nut,t).r));}
void main(){vec2 q=uv;q.x*=resolution.x/resolution.y;mat3 R=rot();vec3 ro=R*vec3(0,0,zoom),rd=R*normalize(vec3(q*0.63,-1.8));vec3 bg=vec3(.0015,.003,.005)*(1.-.2*length(uv));vec3 radiance=vec3(0);vec3 trans=vec3(1);float start=0.;vec3 inv=1./rd;vec3 a=(-vec3(1)-ro)*inv,b=(vec3(1)-ro)*inv;vec3 mn=min(a,b),mx=max(a,b);float near=max(max(mn.x,mn.y),mn.z),far=min(min(mx.x,mx.y),mx.z);
if(far>max(near,0.)){if(mode>0){float t=(slice-ro.z)/rd.z;vec3 p=ro+rd*t;if(t>0.&&all(lessThanEqual(abs(p),vec3(1)))){vec3 f=sampleFields(p);float v=mode==1?(f.x-ds.x)/max(ds.y-ds.x,1e-8):mode==2?(f.y-as.x)/max(as.y-as.x,1e-8):(f.z-ns.x)/max(ns.y-ns.x,1e-8);vec3 c=mode==1?vec3(1,.6,.18):mode==2?vec3(.15,.95,1):vec3(.55,.72,1);bg=mix(vec3(.025,.04,.055),c,pow(clamp(v,0.,1.),.7));}}else{float stepSize=(far-near)/192.;for(int i=0;i<192;i++){vec3 p=ro+rd*(near+(float(i)+.5)*stepSize);if(dot(p/vec3(.86,1.,.78),p/vec3(.86,1.,.78))>1.)continue;vec3 f=sampleFields(p);float d=max(f.x,0.),e=max(f.y,0.);vec3 absorption=d*opacity*vec3(.7,1.6,3.9);vec3 attenuation=exp(-absorption*stepSize*3.);vec3 tp=p*.5+.5;vec3 ag=vec3(texture(act,tp+vec3(voxel.x,0,0)).r-texture(act,tp-vec3(voxel.x,0,0)).r,texture(act,tp+vec3(0,voxel.y,0)).r-texture(act,tp-vec3(0,voxel.y,0)).r,texture(act,tp+vec3(0,0,voxel.z)).r-texture(act,tp-vec3(0,0,voxel.z)).r);float excitation=smoothstep(.08,.65,e)*smoothstep(.025,.20,length(ag/(voxel*40.)));vec3 source=vec3(1.,.53,.12)*d*.42+vec3(.08,.85,1.)*excitation*2.5; radiance+=trans*source*stepSize;trans*=attenuation;if(max(max(trans.r,trans.g),trans.b)<.008)break;}bg=radiance+trans*bg;}}
color=vec4(pow(vec3(1)-exp(-bg*exposure),vec3(.4545)),1);}`;
function shader(type,source){let s=gl.createShader(type);gl.shaderSource(s,source);gl.compileShader(s);if(!gl.getShaderParameter(s,gl.COMPILE_STATUS))throw Error(gl.getShaderInfoLog(s));return s}let program=gl.createProgram();gl.attachShader(program,shader(gl.VERTEX_SHADER,vertex));gl.attachShader(program,shader(gl.FRAGMENT_SHADER,fragment));gl.linkProgram(program);if(!gl.getProgramParameter(program,gl.LINK_STATUS))throw Error(gl.getProgramInfoLog(program));gl.useProgram(program);let buffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,buffer);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array([-1,-1,1,-1,-1,1,-1,1,1,-1,1,1]),gl.STATIC_DRAW);let loc=gl.getAttribLocation(program,'pos');gl.enableVertexAttribArray(loc);gl.vertexAttribPointer(loc,2,gl.FLOAT,false,0,0);const U=name=>gl.getUniformLocation(program,name);let textures=['deposit','activity','nutrient'].map((c,i)=>{let t=gl.createTexture();gl.activeTexture(gl.TEXTURE0+i);gl.bindTexture(gl.TEXTURE_3D,t);gl.texParameteri(gl.TEXTURE_3D,gl.TEXTURE_MIN_FILTER,gl.LINEAR);gl.texParameteri(gl.TEXTURE_3D,gl.TEXTURE_MAG_FILTER,gl.LINEAR);for(let axis of [gl.TEXTURE_WRAP_S,gl.TEXTURE_WRAP_T,gl.TEXTURE_WRAP_R])gl.texParameteri(gl.TEXTURE_3D,axis,gl.CLAMP_TO_EDGE);gl.uniform1i(U(['dep','act','nut'][i]),i);gl.uniform2fv(U(['ds','as','ns'][i]),data.scales[c]);return t});
let angles=[.43,-.18],zoom=3.25,playing=false,last=0,index=10,dirty=true;const decoded={};function load(n){index=n;['deposit','activity','nutrient'].forEach((c,i)=>{let key=c+n;if(!decoded[key]){let b=atob(data.fields[c][n]);decoded[key]=Uint8Array.from(b,x=>x.charCodeAt(0));}gl.activeTexture(gl.TEXTURE0+i);gl.bindTexture(gl.TEXTURE_3D,textures[i]);gl.pixelStorei(gl.UNPACK_ALIGNMENT,1);gl.texImage3D(gl.TEXTURE_3D,0,gl.R8,...data.run.shape,0,gl.RED,gl.UNSIGNED_BYTE,decoded[key]);});let f=data.run.frames[n];$('timestamp').textContent='t '+f.time.toFixed(3);$('solid').textContent=Number(f.solid ?? f.occupied ?? 0).toFixed(1);$('frame').textContent=String(n+1).padStart(2,'0')+' / '+data.run.frames.length;$('time').value=n;dirty=true;}
$('time').max=data.run.frames.length-1;function draw(){let w=Math.round(canvas.clientWidth*Math.min(devicePixelRatio,1.5)),h=Math.round(canvas.clientHeight*Math.min(devicePixelRatio,1.5));if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h;gl.viewport(0,0,w,h);dirty=true;}if(!dirty)return;gl.uniform2f(U('resolution'),w,h);gl.uniform3fv(U('voxel'),data.run.shape.map(n=>1/n));gl.uniform2fv(U('angles'),angles);gl.uniform1f(U('zoom'),zoom);gl.uniform1f(U('opacity'),+$('opacity').value);gl.uniform1f(U('exposure'),+$('exposure').value);gl.uniform1f(U('slice'),+$('slice').value);gl.uniform1i(U('mode'),+$('mode').value);gl.drawArrays(gl.TRIANGLES,0,6);dirty=false;}
$('time').oninput=()=>load(+$('time').value);$('play').onclick=()=>{playing=!playing;$('play').textContent=playing?'Ⅱ Pause':'▶ Play';last=performance.now();};for(let [id,out] of [['opacity','ov'],['exposure','ev'],['slice','sv']]){$(id).oninput=()=>{$(out).textContent=(+$(id).value).toFixed(2);dirty=true;};$(id).oninput();}$('mode').onchange=()=>{dirty=true;let c=['deposit','deposit','activity','nutrient'][+$('mode').value];$('scale').textContent=(data.run.source_channels?.[c]||c)+' range: '+data.scales[c].map(x=>x.toFixed(3)).join(' → ');};$('mode').onchange();let drag=null;canvas.onpointerdown=e=>{drag=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId);};canvas.onpointermove=e=>{if(drag){angles[0]+=(e.clientX-drag[0])*.007;angles[1]+=(e.clientY-drag[1])*.007;drag=[e.clientX,e.clientY];dirty=true;}};canvas.onpointerup=()=>drag=null;canvas.onpointercancel=()=>drag=null;canvas.onwheel=e=>{e.preventDefault();zoom=Math.max(2,Math.min(5,zoom+e.deltaY*.002));dirty=true;};window.onresize=()=>dirty=true;load(index);function tick(t){if(playing&&t-last>550){load((index+1)%data.run.frames.length);last=t;}draw();requestAnimationFrame(tick);}requestAnimationFrame(tick);window.choirViewer={get frame(){return index},get playing(){return playing},scales:data.scales};
}catch(e){$('error').style.display='block';$('error').textContent=e.message;console.error(e);}}
</script></html>'''
title=html_escape.escape(run.get('title','Recorded material'))
if root.name!='choir':
    html=html.replace('Choir — a material remembers',title).replace('Choir <em>/ a material remembers</em>',title)
shape=run['shape']
dimlabel=f"{shape[0]}³" if len(set(shape))==1 else ' × '.join(map(str,shape))
html=html.replace('RECORDED JULIA · CPU · 40³',f"RECORDED JULIA · {run.get('backend','CPU').upper()} · {dimlabel}")
config=run.get('config',{});sim=config.get('simulation',{})
html=html.replace('specimen 017',f"specimen {sim.get('seed',17):03d}")
html=html.replace('seed 17 · dt 0.02',f"seed {sim.get('seed',17)} · dt {sim.get('dt',.02)}")
count=config.get('particles',{}).get('count',0)
html=html.replace('<strong>600</strong>',f'<strong>{count}</strong>')
html=html.replace('21 recorded frames',f"{len(run['frames'])} recorded frames")
if count==0:
    html=html.replace('An excitable medium and recurrent particle policies leave a material trace.','An active medium leaves a material trace.').replace('Fixed recurrent policies; no training or evolution.','Recorded coupled field dynamics; no particle policies in this run.')
html=html.replace('index=10,dirty=true',f'index={min(10,len(run["frames"])-1)},dirty=true')
if not config.get('excitable'):
    html=html.replace('An active medium leaves a material trace. Cyan emphasizes spatial changes in activity inside that accumulated amber. Scrub the recording to watch the balance change.',html_escape.escape(run.get('model','Recorded material fields.'))+' Scrub the recording to inspect its progression.')
    html=html.replace('The simulated domain is a periodic cube; the crop introduces no physical wall.','The crop introduces no physical wall; domain boundaries belong to the recorded model.')
    html=html.replace('Recorded coupled field dynamics; no particle policies in this run.','Recorded field dynamics.')
    if run.get('frames') and 'occupied' in run['frames'][0]:html=html.replace('deposit sum','occupied cells')
if not sim:
    html=html.replace('seed 17 · dt 0.02',html_escape.escape('seed '+str(run.get('random_state',{}).get('seed_u64','unspecified'))))
html=html.replace('x fastest · offline WebGL2','x fastest · offline WebGL2')
if streamed:
    html=html.replace("f.time.toFixed(3)","f.time.toFixed(6)")
    html=html.replace('Active materials / specimen','Recorded material / specimen')
    html=html.replace('Fixed recurrent policies; no training or evolution.','Recorded field dynamics.')
    html=html.replace('Recorded coupled field dynamics; no particle policies in this run.','Recorded field dynamics.')
    html=html.replace('Recorded field dynamics. Playback', 'Recorded field dynamics. Scheduled simulation time is shown. Playback')
    html=html.replace('<strong>0</strong>particles','<strong>'+str(metadata['completed_steps'])+'</strong>steps')
    first=source_channels['deposit'];second=source_channels['activity']
    html=html.replace('Deposit / persistent material',(first or 'Not recorded')+' / amber field')
    html=html.replace('Activity / excitable fronts',(second or 'Not recorded')+' / cyan gradients')
    html=html.replace('AMBER MEMORY · CYAN EXCITATION',html_escape.escape((first or 'none').upper()+' · '+(second or 'none').upper()))
    html=html.replace('Interactive volume of recorded deposit and activity. Drag to rotate.','Interactive recorded volume. Drag to rotate.')
    for role,c in source_channels.items():
        label=html_escape.escape(c or 'not recorded')
        html=html.replace('Slice / '+role,'Slice / '+label)
        if c is None:html=html.replace('value="'+str(['deposit','activity','nutrient'].index(role)+1)+'">Slice','value="'+str(['deposit','activity','nutrient'].index(role)+1)+'" disabled>Slice')
    if first!='deposit':html=html.replace('deposit sum','deposit sum (separate)')
    if metadata['budget']['kind']=='phase':
        html=html.replace('An active medium leaves a material trace. Cyan emphasizes spatial changes in activity inside that accumulated amber. Scrub the recording to watch the balance change.','Amber shows the positive phase; cyan emphasizes gradients in the recorded interface field. Slices show the signed composition and other selected fields directly.')
        html=html.replace(html_escape.escape(run['model'])+' Scrub the recording to inspect its progression.','Amber shows the positive phase; cyan emphasizes gradients in the recorded interface field. Slices show the signed composition and other selected fields directly.')
output=root/'viewer.html';temporary=root/'viewer.html.tmp'
temporary.write_text(html.replace('__DATA__',data));temporary.replace(output);print(output)
