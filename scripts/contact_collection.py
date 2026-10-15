#!/usr/bin/env python3
"""Labeled contact sheets for visual review; never substitutes for source images."""
import argparse,json,pathlib,sys
try: from PIL import Image,ImageDraw
except ImportError:
 sys.path.insert(0,'/tmp/intmat-render-libs'); from PIL import Image,ImageDraw
p=argparse.ArgumentParser();p.add_argument('--root',default='outputs/collection-20261015');p.add_argument('--stage',default='preview');p.add_argument('--family');a=p.parse_args();root=pathlib.Path(a.root)
paths=[]
for selection in root.glob('collection-*-selected.json'):
 data=json.loads(selection.read_text());ids=data if isinstance(data,list) else data['ids']
 for ident in ids:
  path=root/ident/(a.stage+'.png')
  status=path.parent/(a.stage+'-status.json')
  current=status.exists() and json.loads(status.read_text()).get('status')=='complete' and status.stat().st_mtime_ns >= (path.parent/'render-job.json').stat().st_mtime_ns
  if path.exists() and current and (not a.family or json.loads((path.parent/'render-job.json').read_text())['family']==a.family):paths.append(path)
paths=sorted(set(paths))
if not paths:
 print('No current completed '+a.stage+' images'+(' for '+a.family if a.family else ''));sys.exit(0)
cols=6 if len(paths)>12 else 3;w=320;h=350
sheet=Image.new('RGB',(cols*w,((len(paths)+cols-1)//cols)*h),(27,30,35));draw=ImageDraw.Draw(sheet)
for n,path in enumerate(paths):
 im=Image.open(path).convert('RGB');im.thumbnail((w,w));x=n%cols*w;y=n//cols*h;sheet.paste(im,(x+(w-im.width)//2,y));draw.text((x+7,y+w+5),path.parent.name,fill=(235,235,235))
out=root/(a.stage+'-contact'+('-'+a.family if a.family else '')+'.jpg');sheet.save(out,quality=94);print(out.resolve())
