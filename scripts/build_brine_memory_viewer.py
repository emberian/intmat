#!/usr/bin/env python3
"""Adapt four immutable brine fields to the shared compressed volume viewer.

Requires NumPy. Display reduction is nonoverlapping 2×2×2 arithmetic mean,
rounded to Float32, followed by the existing global-per-channel 8-bit mapping.
"""
import argparse,base64,gzip,hashlib,json,pathlib,subprocess,sys,tempfile
from viewer_storage import RECORDING,unpack
import numpy as np

def reduce_field(path):
    v=np.fromfile(path,dtype='<f4');assert v.size==128**3 and np.isfinite(v).all()
    return v.reshape(64,2,64,2,64,2).mean(axis=(1,3,5),dtype=np.float64).astype('<f4').ravel()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run',type=pathlib.Path,nargs='?',default=pathlib.Path('outputs/art/runtime-brine-porous-128-20261016/closed'));p.add_argument('--output',type=pathlib.Path,default=pathlib.Path('site/collection/collection-20261016/brine-left-open/viewer.html'));a=p.parse_args()
    source=json.loads((a.run/'run.json').read_text());assert source['shape']==[128]*3 and len(source['frames'])==4
    frames=[];hashes={}
    with tempfile.TemporaryDirectory(prefix='brine-viewer-') as temporary:
        out=pathlib.Path(temporary)
        for frame in source['frames']:
            name=f'frame-{frame["step"]:05d}';history=next(h for h in source['history'] if h['step']==frame['step'])
            frames.append(dict(name=name,time=frame['time'],solid=history['mineral_mass']))
            for field,alias in [('mineral','deposit'),('phi','activity'),('a','nutrient')]:
                path=a.run/frame['directory']/(field+'.f32');hashes[f'{frame["step"]}/{field}']=hashlib.sha256(path.read_bytes()).hexdigest()
                reduce_field(path).tofile(out/(name+'-'+alias+'.f32'))
        metadata=dict(title=source['title'],shape=[64]*3,native_shape=[128]*3,backend='CPU',frames=frames,config={},
            model='A finite brine leaves a mineral cage inside a changing phase. The mineral starts at zero.',
            source_sha256=source['source_sha256'],capture_sha256=hashes,execution_base_commit=source['execution_base_commit'],
            reduction='Native128³ →64³ by nonoverlapping2×2×2 arithmetic mean, rounded Float32; then global linear8-bit scales per channel. Four actual states; no temporal interpolation.')
        (out/'run.json').write_text(json.dumps(metadata))
        subprocess.run([sys.executable,str(pathlib.Path(__file__).with_name('build_viewer.py')),str(out)],check=True)
        text=(out/'viewer.html').read_text()
    def replace(old,new):
        nonlocal text
        assert old in text,old[:100]
        text=text.replace(old,new)
    replace('<div class="eyebrow">Active materials / specimen 017</div>','<a href="../../#collection-20261016/brine-left-open">← Collection</a><div class="eyebrow">Inside the brine / four saved states</div>')
    replace('RECORDED JULIA · CPU · 64³','NATIVE 128³ · VIEW 64³ · 4 ACTUAL STATES')
    replace('AMBER MEMORY · CYAN EXCITATION','MINERAL · PHASE · REMAINING A')
    replace('Interactive volume of recorded deposit and activity. Drag to rotate.','Interactive mineral, phase and reagent A volume. Drag to rotate.')
    replace('Slice / deposit','Slice / mineral');replace('Slice / activity','Slice / signed phase');replace('Slice / nutrient','Slice / reagent A')
    replace('Volume / amber + cyan','Volume / mineral + phase')
    replace('Slice depth','Cut / slice depth');replace('max="0.95" step="0.025" value="0"','max="1" step="0.025" value="1"')
    replace('if(dot(p/vec3(.86,1.,.78),p/vec3(.86,1.,.78))>1.)continue;','if(p.z>slice)continue;')
    replace('<option value="0">Volume / mineral + phase</option>','<option value="4" selected>Surface / mineral</option><option value="0">Volume / mineral + phase</option>')
    replace('if(mode>0){','if(mode>0 && mode<4){')
    replace("['deposit','deposit','activity','nutrient'][+$('mode').value]","['deposit','deposit','activity','nutrient','deposit'][+$('mode').value]")
    replace('float d=max(f.x,0.),e=max(f.y,0.);', '''float d=max(f.x,0.),e=max(f.y,0.);
if(mode==4){if(d<.30)continue;
float lower=max(near,near+(float(i)-.5)*stepSize),upper=near+(float(i)+.5)*stepSize;
for(int refine=0;refine<8;refine++){float middle=(lower+upper)*.5;vec3 point=ro+rd*middle;
float mineral=mix(ds.x,ds.y,texture(dep,point*.5+.5).r);
if(point.z<=slice && mineral>=.30)upper=middle;else lower=middle;}
bool cutFace=(ro+rd*lower).z>slice;
p=ro+rd*upper;f=sampleFields(p);vec3 tp=p*.5+.5;
vec3 grad=vec3(texture(dep,tp+vec3(voxel.x,0,0)).r-texture(dep,tp-vec3(voxel.x,0,0)).r,texture(dep,tp+vec3(0,voxel.y,0)).r-texture(dep,tp-vec3(0,voxel.y,0)).r,texture(dep,tp+vec3(0,0,voxel.z)).r-texture(dep,tp-vec3(0,0,voxel.z)).r);
vec3 normal=cutFace?vec3(0,0,1):-normalize(grad+vec3(1e-8));if(dot(normal,rd)>0.)normal=-normal;
float light=.22+.68*max(0.,dot(normal,normalize(vec3(-.7,-.8,1.))))+.10*max(0.,dot(normal,normalize(vec3(1.,.4,-.2))));
vec3 pigment=mix(vec3(.69,.30,.13),vec3(.32,.65,.53),clamp((f.z-ns.x)/max(ns.y-ns.x,1e-8),0.,1.));radiance=pigment*light;trans=vec3(0);break;}
''')
    replace('Optical density <span','Density (volume) <span')
    replace('"source_channels":{"deposit":"deposit","activity":"activity","nutrient":"nutrient"}','"source_channels":{"deposit":"mineral","activity":"signed phase","nutrient":"reagent A"}')
    start=text.index('<div class="below">');end=text.index('</footer>',start)+len('</footer>')
    text=text[:start]+'''<div class="below"><div><p>A finite brine leaves an open mineral cage. Turn it, cut through it, or step back through its formation. At the first saved state, mineral is absent.</p><p style="color:#8ba3a2">The surface follows mineral concentration 0.30; its tint follows remaining reagent A. Slices reveal signed phase and the reagent directly. The periodic cube is shown whole, without a viewing crop.</p></div><div><div class="stats"><div><strong id="solid"></strong>native mineral mass</div><div><strong id="frame"></strong>saved state</div><div><strong>128³</strong>native grid</div></div><p style="color:#8ba3a2">The display averages each 2×2×2 block to 64³, then uses fixed 8-bit scales across all four times. Fine pores soften. These browser lights are a separate observation from the photograph.</p><div id="scale" style="font:11px ui-monospace,monospace;color:#8ba3ad"></div></div></div><footer>4 actual states · native 128³ → display 64³ · fixed channel scales · offline WebGL2</footer>'''+text[end:]
    # Lossless differences along x, y and z remove smooth spatial redundancy.
    # Quantization remains the shared viewer's unchanged global8-bit mapping.
    match=RECORDING.search(text);payload=json.loads(match[2]);original=unpack(payload);packed=[]
    for channel,values in original.items():
        for value in values:
            raw=base64.b64decode(value);q=np.frombuffer(raw,dtype=np.uint8).reshape(64,64,64).copy()
            for axis in (2,1,0):
                right=[slice(None)]*3;left=[slice(None)]*3;right[axis]=slice(1,None);left[axis]=slice(None,-1)
                q[tuple(right)]-=q[tuple(left)].copy()
            restored=q.copy()
            for axis in (0,1,2):restored=np.cumsum(restored,axis=axis,dtype=np.uint64).astype(np.uint8)
            assert restored.tobytes()==raw
            packed.append(q.tobytes())
    payload['fields_gzip']=base64.b64encode(gzip.compress(b''.join(packed),compresslevel=9,mtime=0)).decode()
    payload['field_predictor']='axis-difference-u8-v1'
    payload['field_transport']='gzip-u8-channel-frame'
    text=text[:match.start(2)]+json.dumps(payload,separators=(',',':'))+text[match.end(2):]
    replace('data.fields=JSON.parse(await new Response(stream).text());delete data.fields_gzip;',"""const raw=new Uint8Array(await new Response(stream).arrayBuffer());data.fields={};let cursor=0;const count=data.run.shape.reduce((a,b)=>a*b,1);
for(const c of ['deposit','activity','nutrient']){data.fields[c]=[];for(let f=0;f<data.run.frames.length;f++){let value='';for(let k=cursor;k<cursor+count;k+=4096)value+=String.fromCharCode(...raw.subarray(k,Math.min(k+4096,cursor+count)));data.fields[c].push(btoa(value));cursor+=count;}}
if(cursor!==raw.length)throw Error('Recording byte count mismatch');delete data.fields_gzip;""")
    replace('decoded[key]=Uint8Array.from(b,x=>x.charCodeAt(0));',"""decoded[key]=Uint8Array.from(b,x=>x.charCodeAt(0));
if(data.field_predictor==='axis-difference-u8-v1'){const v=decoded[key],nx=data.run.shape[0],ny=data.run.shape[1],plane=nx*ny;
for(let i=0;i<v.length;i++)if(i%nx)v[i]+=v[i-1];
for(let i=0;i<v.length;i++)if(Math.floor(i/nx)%ny)v[i]+=v[i-nx];
for(let i=plane;i<v.length;i++)v[i]+=v[i-plane];}""")
    assert len(text.encode())<1_000_000,f'Viewer exceeds1MB: {len(text.encode())}'
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(text);print(a.output,len(text.encode()),'bytes')
if __name__=='__main__':main()
