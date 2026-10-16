"""CPU photograph actual stored fracture states, without temporal interpolation."""
import argparse,json,runpy,sys,time,hashlib,shutil
import bpy
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('capture',type=Path);p.add_argument('output',type=Path);p.add_argument('--prototype',action='store_true');p.add_argument('--reuse',type=Path);p.add_argument('--size',type=int,default=640);p.add_argument('--samples',type=int,default=16)
a=p.parse_args(sys.argv[sys.argv.index('--')+1:]);files=sorted((a.capture/'states').glob('*.json'));assert len(files)==40
if a.prototype:files=[files[0],files[20],files[-1]]
a.output.mkdir(parents=True,exist_ok=True);reuse={};renderer=Path(__file__).with_name('render_fault_remembers.py')
if a.reuse:
 meta=json.loads((a.reuse/'frames.json').read_text())
 assert (meta['size'],meta['samples'],meta['renderer_sha256'])==(a.size,a.samples,hashlib.sha256(renderer.read_bytes()).hexdigest()), 'Reused photography settings differ'
 for rec in meta['frames']:reuse[rec['source_sha256']]=(rec,a.reuse/rec['image'])
records=[]
for i,path in enumerate(files):
 raw=json.loads(path.read_text());start=time.perf_counter();target=a.output/f'{i:04d}.png'
 sys.argv=[str(renderer),'--',str(path),str(target),'--size',str(a.size),'--samples',str(a.samples),'--dose-max','.04419005268498425','--ortho-scale','3.2']
 source_hash=hashlib.sha256(path.read_bytes()).hexdigest()
 if source_hash in reuse:
  rec,image=reuse[source_hash]
  if image.resolve()!=target.resolve():shutil.copy2(image,target);shutil.copy2(image.with_suffix('.json'),target.with_suffix('.json'))
 else:runpy.run_path(str(renderer),run_name='__main__')
 bpy.data.orphans_purge(do_recursive=True)
 record={'index':i,'step':raw['step'],'time':raw['time'],'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'image':target.name,'elapsed_seconds':time.perf_counter()-start};records.append(record)
 (a.output/'frames.json').write_text(json.dumps({'actual_frames':True,'temporal_interpolation':False,'dose_optical_max':.04419005268498425,'renderer_sha256':hashlib.sha256(renderer.read_bytes()).hexdigest(),'size':a.size,'samples':a.samples,'frames':records},indent=2))
 print('COMPLETE',json.dumps(record),flush=True)
