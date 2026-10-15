#!/usr/bin/env python3
"""Bounded read-only observer of SHA-verified completed native-garden frames."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import tempfile
import time

ROOT=Path(__file__).resolve().parents[1]

def atomic_json(path,value):
    temporary=path.with_suffix(path.suffix+'.pending')
    temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)

def completed_entries(source,manifest,known=()):
    entries=[]
    for entry in manifest.get('frames',[]):
        if entry.get('sha256') in known:continue
        relative=PurePosixPath(entry.get('file',''))
        if relative.is_absolute() or '..' in relative.parts or not relative.parts or relative.parts[0]!='frames' or relative.suffix!='.jls':continue
        path=source.joinpath(*relative.parts)
        if not path.is_file() or not path.resolve().is_relative_to(source.resolve()):continue
        if path.stat().st_size!=entry.get('bytes'):continue
        if hashlib.sha256(path.read_bytes()).hexdigest()!=entry.get('sha256'):continue
        entries.append((entry,path))
    return entries

class Observer:
    def __init__(self,source,destination,budget=256_000_000,julia='julia',export_helper=None):
        self.source=Path(source).resolve();self.destination=Path(destination).resolve();self.budget=budget;self.julia=julia;self.export_helper=Path(export_helper) if export_helper else ROOT/'scripts/export_native_observatory_frame.jl'
        self.destination.mkdir(parents=True,exist_ok=True)
        self.index={'schema':'native-garden-observatory-1','frames':[],'observer':{},'source_name':self.source.name}
        existing=self.destination/'observatory.json'
        if existing.exists():
            self.index=json.loads(existing.read_text())
            if self.index.get('source_name')!=self.source.name:raise ValueError('destination belongs to another run')
        identity=self.index.get('source_identity')
        if identity and identity.get('resolved_path')!=str(self.source):
            raise ValueError('destination belongs to another source path')
        # Legacy indexes can be bound only when every preserved observation is
        # explicitly recorded by the supplied source, never from its basename alone.
        if not identity and self.index['frames']:
            try:
                manifest=json.loads((self.source/'run.json').read_text())
                recorded={f['sha256'] for f in manifest['frames']}
                if not all(f['source_sha256'] in recorded for f in self.index['frames']):
                    raise ValueError('legacy observations do not belong to supplied source')
            except (OSError,KeyError,ValueError) as error:
                raise ValueError('cannot establish legacy observer source identity') from error
        self.index.setdefault('source_identity',{'resolved_path':str(self.source)})
        previous=self.index['observer'].pop('stop_reason',None)
        if previous:self.index['observer']['previous_session_stop_reason']=previous
        self.index['observer'].pop('error',None)
        shutil.copyfile(ROOT/'web/native_observatory.html',self.destination/'index.html')
        self.last_manifest=None
    def size(self):return sum(p.stat().st_size for p in self.destination.rglob('*') if p.is_file())
    def scan(self):
        try:
            manifest=json.loads((self.source/'run.json').read_text())
            if not isinstance(manifest.get('frames'),list):raise ValueError('no completed frame entries')
            self.last_manifest=manifest
        except (OSError,ValueError):
            self.index['observer'].update(status='waiting for a complete manifest; retained last good observations',updated_unix=time.time())
            atomic_json(self.destination/'observatory.json',self.index);return False
        identity=self.index['source_identity']
        first=manifest['frames'][0].get('sha256') if manifest['frames'] else None
        if identity.get('first_frame_sha256') and first!=identity['first_frame_sha256']:
            raise ValueError('source run first-frame identity changed')
        if first:identity['first_frame_sha256']=first
        known={f['source_sha256'] for f in self.index['frames']}
        for entry,path in completed_entries(self.source,manifest,known):
            if entry['sha256'] in known:continue
            if self.size()+entry['bytes']*2+8192>self.budget:
                self.index['observer']['stop_reason']='observer artifact budget';break
            name=f"frame-{entry['step']:09d}-{entry['side']}-{entry['sha256'][:10]}"
            temporary=Path(tempfile.mkdtemp(prefix='.export-',dir=self.destination));target=temporary/'payload'
            try:
                subprocess.run([self.julia,'--startup-file=no','--project='+str(ROOT),str(self.export_helper),str(path),entry['sha256'],str(target)],check=True,timeout=120)
                metadata=json.loads((target/'frame.json').read_text())
                if metadata['source_sha256']!=entry['sha256'] or metadata['step']!=entry['step'] or metadata['shape']!=[entry['side']]*3 or abs(metadata['time']-entry['time'])>1e-9:
                    raise ValueError('export metadata disagrees with completed manifest entry')
            except (OSError,ValueError,subprocess.SubprocessError) as error:
                self.index['observer'].update(stop_reason='export failed; retained last good observations',error=str(error))
                break
            # Only complete verified exports become visible in the observer manifest.
            if self.size()>self.budget:
                self.index['observer']['stop_reason']='observer artifact budget';break
            final=self.destination/name;target.replace(final);temporary.rmdir()
            self.index['frames'].append({'file':name+'/frame.json','time':metadata['time'],'step':metadata['step'],'side':metadata['shape'][0],
                'native_shape':metadata['native_shape'],'population':metadata['population'],'source_sha256':entry['sha256']})
            known.add(entry['sha256'])
        self.index['frames'].sort(key=lambda f:(f['time'],f['side']))
        self.index.update(run_complete=bool(manifest.get('complete')),deadline_unix=manifest.get('deadline_unix'),
            source_commit=manifest.get('source_commit'),limits=manifest.get('limits'),native_shape=manifest.get('config',{}).get('grid',{}).get('shape'),
            physical_dt=manifest.get('config',{}).get('native',{}).get('dt'))
        self.index['observer'].update(status='completed frames only',updated_unix=time.time(),artifact_budget_bytes=self.budget,stored_bytes=self.size())
        atomic_json(self.destination/'observatory.json',self.index);return True

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run',type=Path);p.add_argument('destination',type=Path)
    p.add_argument('--follow',action='store_true');p.add_argument('--max-wall-seconds',type=float,default=600);p.add_argument('--poll-seconds',type=float,default=5)
    p.add_argument('--artifact-budget-bytes',type=int,default=256_000_000);p.add_argument('--julia',default='julia');p.add_argument('--export-helper',type=Path);a=p.parse_args()
    if a.max_wall_seconds<=0 or a.poll_seconds<=0 or a.artifact_budget_bytes<=0:p.error('positive bounds required')
    observer=Observer(a.run,a.destination,a.artifact_budget_bytes,a.julia,a.export_helper);wall_deadline=time.time()+a.max_wall_seconds
    while True:
        observer.scan();deadline=min(wall_deadline,observer.index.get('deadline_unix') or wall_deadline)
        if not a.follow or observer.index.get('run_complete') or observer.index['observer'].get('stop_reason') or time.time()>=deadline:break
        time.sleep(min(a.poll_seconds,max(0,deadline-time.time())))
    observer.index['observer'].setdefault('stop_reason','run complete' if observer.index.get('run_complete') else 'single scan' if not a.follow else 'observer/run wall deadline')
    atomic_json(a.destination/'observatory.json',observer.index)
    print(json.dumps({'frames':len(observer.index['frames']),'destination':str(a.destination),'status':observer.index['observer']}))
if __name__=='__main__':main()
