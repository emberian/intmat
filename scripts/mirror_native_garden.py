#!/usr/bin/env python3
"""Read-only bounded SSH mirror of completed, SHA-verified garden frames."""
import argparse,hashlib,json,shlex,subprocess,time
from pathlib import Path,PurePosixPath
from observe_native_garden import atomic_json

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('host');p.add_argument('remote');p.add_argument('destination',type=Path)
    p.add_argument('--deadline-unix',type=float,required=True);p.add_argument('--budget',type=int,default=1_000_000_000);p.add_argument('--poll',type=float,default=45)
    a=p.parse_args();a.destination.mkdir(parents=True,exist_ok=True)
    identity={'host':a.host,'remote':a.remote};ip=a.destination/'mirror-identity.json'
    if ip.exists() and json.loads(ip.read_text())!=identity:raise ValueError('mirror destination belongs to another source')
    atomic_json(ip,identity);status={'source':identity,'budget_bytes':a.budget,'deadline_unix':a.deadline_unix}
    def remote(relative):return ['ssh','-o','BatchMode=yes','-o','ConnectTimeout=15',a.host,'cat -- '+shlex.quote(a.remote+'/'+relative)]
    while time.time()<a.deadline_unix:
        try:
            raw=subprocess.check_output(remote('run.json'),timeout=30);manifest=json.loads(raw)
            if not isinstance(manifest['frames'],list):raise ValueError('invalid frame manifest')
            complete=[];budget_stop=False
            for entry in manifest['frames']:
                rel=PurePosixPath(entry['file'])
                if rel.is_absolute() or '..' in rel.parts or rel.parts[0]!='frames' or rel.suffix!='.jls':raise ValueError('unsafe frame path')
                target=a.destination.joinpath(*rel.parts)
                if target.exists():
                    if target.stat().st_size!=entry['bytes'] or hashlib.sha256(target.read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('retained frame changed')
                else:
                    used=sum(f.stat().st_size for f in a.destination.rglob('*') if f.is_file())
                    if used+entry['bytes']+len(raw)>a.budget:budget_stop=True;break
                    target.parent.mkdir(parents=True,exist_ok=True);partial=target.with_suffix('.partial')
                    with partial.open('wb') as out:subprocess.run(remote(str(rel)),stdout=out,check=True,timeout=120)
                    if partial.stat().st_size!=entry['bytes'] or hashlib.sha256(partial.read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('downloaded frame hash mismatch')
                    partial.replace(target)
                complete.append(entry)
            evidence=a.destination/'manifests'/('run-'+hashlib.sha256(raw).hexdigest()+'.json');evidence.parent.mkdir(exist_ok=True)
            if not evidence.exists():evidence.write_bytes(raw)
            manifest['frames']=complete;manifest['complete']=bool(manifest.get('complete')) and len(complete)==len(json.loads(raw)['frames'])
            atomic_json(a.destination/'run.json',manifest)
            status.update(frames=len(complete),updated_unix=time.time(),status='completed verified frames only')
            if budget_stop:status['stop_reason']='mirror artifact budget'
            elif manifest['complete']:status['stop_reason']='remote run complete'
            atomic_json(a.destination/'mirror-status.json',status)
            if status.get('stop_reason'):break
        except (OSError,ValueError,KeyError,subprocess.SubprocessError) as error:
            status.update(error=str(error),updated_unix=time.time(),status='retaining last good mirror');atomic_json(a.destination/'mirror-status.json',status)
        time.sleep(min(a.poll,max(0,a.deadline_unix-time.time())))
    status.setdefault('stop_reason','hard deadline');atomic_json(a.destination/'mirror-status.json',status)
if __name__=='__main__':main()
