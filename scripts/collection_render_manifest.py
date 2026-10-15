#!/usr/bin/env python3
"""Verify selected collection photographs and emit their portable review manifest."""
import argparse,hashlib,json,pathlib,sys
p=argparse.ArgumentParser();p.add_argument('--root',default='outputs/collection-20261015');p.add_argument('--require-complete',action='store_true');a=p.parse_args();root=pathlib.Path(a.root).resolve();ids=[]
for selection in sorted(root.glob('collection-*-selected.json')):
 data=json.loads(selection.read_text());ids.extend(data if isinstance(data,list) else data['ids'])
records=[]
for ident in sorted(set(ids)):
 d=root/ident;record={'id':ident,'status':'pending'};job=d/'render-job.json';status=d/'final-status.json';approval=d/'art-reviewed.json';image=d/'final.png'
 if job.exists():
  j=json.loads(job.read_text());record.update(family=j['family'],title=j['title'])
 if status.exists():
  s=json.loads(status.read_text());record.update(status=s['status'],render_seconds=s.get('render_seconds'))
  if s['status']=='complete' and image.exists():
   raw=image.read_bytes();sha=hashlib.sha256(raw).hexdigest();width=int.from_bytes(raw[16:20],'big');height=int.from_bytes(raw[20:24],'big');optical=json.loads((d/'final.png.optical.json').read_text());review=json.loads(approval.read_text()) if approval.exists() else {}
   valid=(sha==s.get('sha256')==review.get('sha256') and review.get('status')=='approved' and [width,height]==[1600,1600] and optical['settings']['samples']==96 and status.stat().st_mtime_ns>=job.stat().st_mtime_ns)
   record.update(status='approved' if valid else 'awaiting_review',image=str(image),sha256=sha,dimensions=[width,height],samples=optical['settings']['samples'],adaptive_sampling=optical.get('adaptive_sampling'),device=optical.get('render_device'),scene=str(d/'final.blend'),vdb=str(d/'final.vdb'),log=str(d/'final.log'),optical=str(d/'final.png.optical.json'))
 records.append(record)
manifest={'selected_count':len(records),'approved_count':sum(r['status']=='approved' for r in records),'records':records}
out=root/'render-review-manifest.json';out.write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps({'manifest':str(out),'selected':len(records),'approved':manifest['approved_count']}))
if a.require_complete and (len(records)!=36 or manifest['approved_count']!=36):sys.exit(1)
