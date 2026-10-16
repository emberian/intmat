#!/usr/bin/env python3
"""Extend the recorded-field/ray viewer with actual absorbed-power captures."""
import argparse,base64,gzip,hashlib,json,pathlib,re,subprocess,sys
import numpy as np
from viewer_storage import RECORDING,unpack
p=argparse.ArgumentParser();p.add_argument('run',type=pathlib.Path);p.add_argument('--tag',default='feedback-final');a=p.parse_args();root=a.run
subprocess.run([sys.executable,str(pathlib.Path(__file__).with_name('build_refractive_viewer.py')),str(root),'--tag',a.tag],check=True)
s=(root/'viewer.html').read_text();m=RECORDING.search(s);data=json.loads(m[2]);fields=unpack(data);run=json.loads((root/'run.json').read_text());n=run['shape'][0];factor=n//64;captures=[];hashes=[]
for frame in run['frames']:
    path=root/frame['directory']/'dose.f32';dose=np.fromfile(path,dtype='<f4').reshape((n,n,n),order='F');reduced=dose.reshape(64,factor,64,factor,64,factor).mean(axis=(1,3,5),dtype=np.float64)
    captures.append(np.log1p(reduced));hashes.append(hashlib.sha256(path.read_bytes()).hexdigest())
hi=max(float(d.max()) for d in captures);data['scales']['activity']=[0,hi]
for frame,raw in zip(data['run']['frames'],run['frames']):frame['dose_time']=None if raw['step']==0 else ((raw['step']-1)//run['params']['refresh_steps'])*run['params']['refresh_steps']*float(run['dt'])
fields['activity']=[base64.b64encode(np.rint(d/hi*255).astype('u1').tobytes(order='F')).decode() for d in captures]
data['fields_gzip']=base64.b64encode(gzip.compress(json.dumps(fields,separators=(',',':')).encode(),mtime=0)).decode();data.pop('fields',None)
data['run']['source_channels']['activity']='log1p absorbed power density';data['run']['dose_sha256']=hashes;data['run']['energy']=run['energy'];data['run']['dose_interpretation']='Block mean of native absorbed-power density, then log1p and fixed 8-bit scaling. State0 precedes illumination; later captures retain power density used in the preceding interval. Material volume units differ from optical length.'
s=s[:m.start(2)]+json.dumps(data,separators=(',',':'))+s[m.end(2):]
s=s.replace('../../#collection-20261016/refractive-memory','../../#collection-20261016/illuminated-memory').replace('A mineral changes the path of light','Light changes the mineral that bends it').replace('AMBER MINERAL · CYAN TRACED PATHS','AMBER MINERAL · CYAN ABSORPTION / PATHS').replace('Paths / fixed launch apertures','Absorbed power / logarithmic display')
s=s.replace('<div class="below">','<label>Recorded fields<select id="fieldlayer"><option value="0">Solid + absorbed power</option><option value="1">Solid only</option><option value="2">Absorbed power only</option></select></label><p id="doseclock"></p><div class="below">')
s=s.replace('uniform int mode;','uniform int mode,fieldLayer;').replace('float d=max(f.x,0.),e=max(f.y,0.);','float d=fieldLayer==2?0.:max(f.x,0.),e=fieldLayer==1?0.:max(f.y,0.);')
s=s.replace("gl.uniform1i(U('mode'),+$('mode').value);","gl.uniform1i(U('mode'),+$('mode').value);gl.uniform1i(U('fieldLayer'),+$('fieldlayer').value);")
s=s.replace("loadRays(n);dirty=true;}","loadRays(n);$('doseclock').textContent=f.dose_time===null?'Initial state: no illumination has been applied; absorbed power is zero.':'Absorbed-power field held from optical refresh t '+f.dose_time.toFixed(3)+'. Geometric path overlay is retraced at the displayed material time.';dirty=true;}")
s=s.replace("for(let id of ['raymode','raystride','cliprays'])","for(let id of ['raymode','raystride','cliprays','fieldlayer'])")
s=s.replace('<option value="1">Slice / solid</option>','<option value="1">Slice / solid</option><option value="2">Slice / log absorbed power</option>')
s=s.replace('Five recorded mineral states, one fixed set of launch apertures. Compare bent paths with their constant-index reference; turn the material or open its interior.','Five recorded states of a light-remodeled mineral. Absorbed power drives its recurrent interior, which changes its next illumination. Cyan field accents reveal the recorded absorption structure; thin paths remain geometric ray observations.')
s=s.replace('The photograph uses native mineral contours.','The photograph uses native mineral and absorbed-power contours. Cyan field accents here show gradients of log1p absorbed power density, averaged to64³ then fixed-scale8-bit encoded; they are a display of recorded energy deposition, not luminescence. The captured power field is held from the preceding optical refresh, while paths are independently retraced at the captured state. The photograph uses native mineral contours.')
energy=run['energy'];sentence=f"Incident energy {energy['incident']:,.2f}; absorbed {energy['absorbed']:,.2f}; transmitted {energy['emergent']:,.2f}; unresolved {energy['unresolved']:,.2f}, in declared arbitrary units. These energy totals cover the whole coupled run; the receiver image is a separate equal-weight geometric probe at this snapshot."
s=s.replace('<h2 style="font:30px Georgia;color:#f6ead4">Where the rays arrive</h2>','<h2 style="font:30px Georgia;color:#f6ead4">Where the rays arrive</h2><p>'+sentence+'</p>')
s=s.replace('</footer>','</footer><p style="color:#91a5ad">Fine outgoing paths are sensitive to the integration step; unresolved paths remain in the accounting. The material response and absorbed-energy balances have separate controls. <a href="https://github.com/emberian/intmat/blob/dev/docs/illuminated-memory.md">Process and numerical limits</a>.</p>')
(root/'viewer.html').write_text(s);print('ILLUMINATED_VIEWER',len(s.encode()))
