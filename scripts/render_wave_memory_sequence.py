"""Render actual stored Julia states with one fixed camera and optical transfer."""
import argparse,hashlib,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('capture',type=Path);p.add_argument('output',type=Path);p.add_argument('--indices');p.add_argument('--size',type=int,default=640);p.add_argument('--samples',type=int,default=24);a=p.parse_args()
root=Path(__file__).resolve().parents[1];renderer=root/'scripts/render_wave_memory_frame.py';pin=hashlib.sha256(renderer.read_bytes()).hexdigest();run=json.loads((a.capture/'run.json').read_text());states=run['frames'][1:];assert len(states)==40
indices=list(map(int,a.indices.split(','))) if a.indices else list(range(40));a.output.mkdir(parents=True,exist_ok=True);path=a.output/'frames.json';record=json.loads(path.read_text()) if path.exists() else {'renderer_sha256':pin,'size':a.size,'samples':a.samples,'camera':'fixed','frames':{}}
assert (record['renderer_sha256'],record['size'],record['samples'])==(pin,a.size,a.samples)
for i in indices:
 state=states[i];volume=a.capture/state['directory'];image=a.output/f'{i:05d}.png';hashes={q.name:hashlib.sha256(q.read_bytes()).hexdigest() for q in sorted(volume.iterdir()) if q.is_file()}
 if str(i) in record['frames']:
  assert record['frames'][str(i)]['source_fields']==hashes and image.is_file();continue
 subprocess.run(['/Applications/Blender.app/Contents/MacOS/Blender','-b','--python',str(renderer),'--',str(volume),str(image),'--size',str(a.size),'--samples',str(a.samples),'--device','METAL','--threads','2','--no-floor'],check=True)
 record['frames'][str(i)]={'image':image.name,'step':state['step'],'time':state['time'],'source_fields':hashes,'image_sha256':hashlib.sha256(image.read_bytes()).hexdigest()};path.write_text(json.dumps(record,indent=2));print('RECORDED',i,state['step'],flush=True)
