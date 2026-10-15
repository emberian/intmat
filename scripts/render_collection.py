#!/usr/bin/env python3
"""Resumable bounded Blender queue. Each photograph retains its own VDB, scene and log."""
import argparse, concurrent.futures, hashlib, json, os, pathlib, subprocess, time
P=argparse.ArgumentParser(); P.add_argument('--root',default='outputs/collection-20261015'); P.add_argument('--stage',choices=['preview','final'],default='preview'); P.add_argument('--jobs',nargs='*'); P.add_argument('--workers',type=int,default=2); P.add_argument('--threads',type=int,default=4); P.add_argument('--watch',action='store_true'); P.add_argument('--device',choices=['CPU','METAL'],default='CPU'); P.add_argument('--selected',action='store_true'); P.add_argument('--approved',action='store_true'); A=P.parse_args()
ROOT=pathlib.Path(A.root).resolve(); REPO=pathlib.Path(__file__).resolve().parent.parent
BLENDER='/Applications/Blender.app/Contents/MacOS/Blender'
def run(path):
 j=json.loads(path.read_text()); d=path.parent; out=d/(A.stage+'.png'); status=d/(A.stage+'-status.json')
 args=j.get('args',[])[:]
 for flag in ['--bond-graph','--population']:
  if flag in args:
   i=args.index(flag)+1; args[i]=str((d/args[i]).resolve())
 command=[BLENDER,'-b','-t',str(A.threads),'--python',str(REPO/'scripts/render_blender.py'),'--',str(d/j.get('volume','volume')),str(out),'--device',A.device,'--size',str(400 if A.stage=='preview' else j.get('size',1600)),'--samples',str(16 if A.stage=='preview' else j.get('samples',96))]+args
 volume_dir=d/j.get('volume','volume')
 source_signature=[(f.name,f.stat().st_size,f.stat().st_mtime_ns) for f in sorted(volume_dir.glob('*')) if f.suffix in ('.f32','.json') and f.name!='optical-mapping.json']
 fingerprint=hashlib.sha256(json.dumps([command,source_signature]).encode()+path.read_bytes()+(REPO/'scripts/render_blender.py').read_bytes()).hexdigest()
 if status.exists() and out.exists():
  old=json.loads(status.read_text())
  if old.get('status')=='complete' and old.get('fingerprint')==fingerprint:return old
 start=time.time(); record={'id':j['id'],'family':j['family'],'stage':A.stage,'image':str(out),'status':'running','fingerprint':fingerprint,'command':command}
 status.write_text(json.dumps(record,indent=2)); env=os.environ.copy(); env.update(OMP_NUM_THREADS=str(A.threads),OPENBLAS_NUM_THREADS='1')
 with (d/(A.stage+'.log')).open('w') as log: r=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env)
 record.update(status='complete' if r.returncode==0 and out.exists() else 'failed',returncode=r.returncode,render_seconds=round(time.time()-start,3))
 if out.exists():
  record['sha256']=hashlib.sha256(out.read_bytes()).hexdigest(); record['dimensions']=[int.from_bytes(out.read_bytes()[16:20],'big'),int.from_bytes(out.read_bytes()[20:24],'big')]
 status.write_text(json.dumps(record,indent=2)); print(json.dumps(record),flush=True); return record
seen=set()
with concurrent.futures.ThreadPoolExecutor(max_workers=min(A.workers,2)) as pool:
 while True:
  selected=set()
  for manifest in ROOT.glob('collection-*-selected.json'):
   data=json.loads(manifest.read_text()); selected.update(data if isinstance(data,list) else data.get('ids',[]))
  paths=[p for p in sorted(ROOT.glob('*/render-job.json')) if (not A.jobs or p.parent.name in A.jobs) and (not A.selected or p.parent.name in selected) and (not A.approved or (p.parent/'art-approved.json').exists()) and (p.parent/'volume/volume.json').exists() and not p.parent.name.startswith('pilot-')]
  pending=[p for p in paths if (str(p),p.stat().st_mtime_ns) not in seen]
  for p in pending: seen.add((str(p),p.stat().st_mtime_ns))
  list(pool.map(run,pending))
  records=[]
  for p in ROOT.glob('*/'+A.stage+'-status.json'):
   records.append(json.loads(p.read_text()))
  (ROOT/(A.stage+'-manifest.json')).write_text(json.dumps(records,indent=2))
  if not A.watch:break
  time.sleep(10)
