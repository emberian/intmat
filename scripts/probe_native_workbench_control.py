#!/usr/bin/env python3
"""Exercise the bounded file protocol of an already launched fresh Julia worker."""
import hashlib,json,pathlib,sys,time
root=pathlib.Path(sys.argv[1]);report=pathlib.Path(sys.argv[2]);deadline=time.monotonic()+450

def load(name):
 try:return json.loads((root/name).read_text())
 except (OSError,ValueError):return None

def wait(predicate):
 while time.monotonic()<deadline:
  value=predicate()
  if value:return value
  time.sleep(.05)
 raise TimeoutError('bounded worker protocol deadline')

manifest=wait(lambda:load('worker.json'))
assert manifest['status']=='ready'
assert all(b['step']==0 for b in manifest['branches'].values())
rows=[]
commands=[('describe',{}),('step',{'branch':'both','steps':16}),
 ('checkpoint',{'branch':'both'}),('pause',{}),('observe',{'branch':'off'}),
 ('step',{'branch':'sensing','steps':8}),('checkpoint',{'branch':'sensing'}),
 ('step',{'branch':'both','steps':257}),('stop',{})]
for number,(action,values) in enumerate(commands,1):
 command=dict(id=f'protocol-{number}',action=action,requested_at=time.time(),**values)
 temporary=root/'command.json.pending';temporary.write_text(json.dumps(command));temporary.replace(root/'command.json')
 start=time.monotonic()
 ack=wait(lambda: (a if a and a.get('id')==command['id'] else None) if (a:=load('control.json')) else None)
 assert ack['ok']==(number!=8),ack
 rows.append(dict(command=command,ack=ack,ack_seconds=time.monotonic()-start))
end=wait(lambda: (m if m and m.get('status')=='stopped' else None) if (m:=load('worker.json')) else None)
assert end['branches']['off']['step']==16 and end['branches']['sensing']['step']==24
frames={}
for branch in ('off','sensing'):
 index=json.loads((root/'branches'/branch/'observatory.json').read_text())
 frames[branch]=[]
 for entry in index['frames']:
  path=root/'branches'/branch/entry['file'];meta=json.loads(path.read_text());binary=path.parent/meta['binary']
  assert hashlib.sha256(binary.read_bytes()).hexdigest()==meta['binary_sha256']
  frames[branch].append(dict(time=meta['time'],step=meta['step'],population=meta['population'],shape=meta['shape'],binary_bytes=binary.stat().st_size,binary_sha256=meta['binary_sha256']))
report.write_text(json.dumps(dict(initial=manifest,final=end,commands=rows,frames=frames),indent=2)+'\n')
print(json.dumps(dict(status=end['status'],steps={k:v['step'] for k,v in end['branches'].items()},commands=len(rows),frames={k:len(v) for k,v in frames.items()})))
