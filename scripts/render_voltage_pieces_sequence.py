"""Fixed-camera photographs of actual electrically consistent fracture observations."""
import argparse,hashlib,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('capture',type=Path);p.add_argument('output',type=Path);p.add_argument('--indices');p.add_argument('--size',type=int,default=640);p.add_argument('--samples',type=int,default=24);a=p.parse_args()
root=Path(__file__).resolve().parents[1];renderer=root/'scripts/render_voltage_pieces.py';pin=hashlib.sha256(renderer.read_bytes()).hexdigest();capture=json.loads((a.capture/'capture.json').read_text());states=capture['frames'];indices=list(map(int,a.indices.split(','))) if a.indices else range(len(states));a.output.mkdir(parents=True,exist_ok=True);path=a.output/'frames.json';dose=capture['global_dose_optical_max']
record=json.loads(path.read_text()) if path.exists() else {'renderer_sha256':pin,'base_sha256':capture['base_sha256'],'size':a.size,'samples':a.samples,'global_dose_optical_max':dose,'camera':'fixed','frames':{}}
assert (record['renderer_sha256'],record['size'],record['samples'],record['global_dose_optical_max'])==(pin,a.size,a.samples,dose)
assert hashlib.sha256((a.capture/'base.json').read_bytes()).hexdigest()==record['base_sha256']
for i in indices:
 s=states[i];state=a.capture/s['file'];assert hashlib.sha256(state.read_bytes()).hexdigest()==s['sha256'];image=a.output/f'{i:05d}.png'
 if str(i) in record['frames']:assert image.is_file() and record['frames'][str(i)]['source_sha256']==s['sha256'];continue
 subprocess.run(['/Applications/Blender.app/Contents/MacOS/Blender','-b','--python',str(renderer),'--',str(state),str(image),'--base',str(a.capture/'base.json'),'--dose-max',str(dose),'--size',str(a.size),'--samples',str(a.samples),'--device','METAL'],check=True)
 record['frames'][str(i)]={'image':image.name,'step':s['step'],'time':s['time'],'source_sha256':s['sha256'],'image_sha256':hashlib.sha256(image.read_bytes()).hexdigest()};path.write_text(json.dumps(record,indent=2));print('RECORDED',i,s['step'],flush=True)
