#!/usr/bin/env python3
"""Plot saved native-field slices with a shared exposure; no model execution."""
import json, math, sys, hashlib
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
src,out=map(Path,sys.argv[1:3]);assert not out.exists()
x=json.loads(src.read_text())
font='/System/Library/Fonts/Supplemental/Arial.ttf'
def f(n):return ImageFont.truetype(font,n)
im=Image.new('RGB',(1280,1060),'#fbfaf6');d=ImageDraw.Draw(im)
d.text((36,25),'A paid patch at t60: actual spatial fields',font=f(28),fill='#283238')
d.text((36,65),'Two central z planes averaged · both chemical species summed · same color exposure across all arms',font=f(16),fill='#283238')
keys=['immobile_full','immobile_catalytic','mobile_full','mobile_catalytic']
labels=['Immobile raw / matrix','Immobile raw / catalyst only','Mobile raw / matrix','Mobile raw / catalyst only']
fields=['raw','free','bound'];caps=[6,250,8]
colors=[(84,111,164),(196,94,48),(82,137,108)]
for col,(field,cap,color) in enumerate(zip(fields,caps,colors)):
 x0=320+col*310
 d.text((x0,104),f'{field.capitalize()} concentration',font=f(19),fill='#283238')
 d.text((x0,130),f'log1p scale, 0 to {cap:g}',font=f(14),fill='#45535a')
 for row,key in enumerate(keys):
  y0=166+row*205;a=x[key][field];n=len(a)
  heat=Image.new('RGB',(n,n));pixels=heat.load()
  for y in range(n):
   for xx in range(n):
    v=min(1,max(0,math.log1p(a[y][xx])/math.log1p(cap)))
    pixels[xx,n-1-y]=tuple(round(249*(1-v)+c*v) for c in color)
  im.paste(heat.resize((184,184),Image.Resampling.NEAREST),(x0,y0))
  # Physical location of the fixed builder and near/far sensor centers.
  cx=x0+92;cy=y0+92
  d.ellipse((cx-16.56,cy-16.56,cx+16.56,cy+16.56),outline='#18242b',width=1)
  d.line((cx-4,cy,cx+4,cy),fill='#18242b',width=1);d.line((cx,cy-4,cx,cy+4),fill='#18242b',width=1)
  far=cx+36.8;d.rectangle((far-2,cy-2,far+2,cy+2),fill='#18242b')
  if col==0:
   d.text((36,y0+59),labels[row],font=f(18),fill='#283238')
   d.text((36,y0+87),'2 × 2 physical extent',font=f(14),fill='#45535a')
 # Vertical legend
 for yy in range(184):
  t=1-yy/183;d.line((x0+198,166+yy,x0+212,166+yy),fill=tuple(round(249*(1-t)+c*t) for c in color))
 d.text((x0+220,166),str(cap),font=f(12),fill='#283238');d.text((x0+220,337),'0',font=f(12),fill='#283238')
d.text((36,1000),'Circle: exchange support R=0.18. Cross: builder sensor. Square: sensor 0.4 away. Frozen builder does not consume food.',font=f(15),fill='#283238')
d.text((36,1026),'Matrix-off retains body porosity. Colors show measured voxel values; logarithmic exposure is explicit and fixed per field.',font=f(15),fill='#283238')
im.save(out)
out.with_suffix('.identity.json').write_text(json.dumps({'slice_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'plot_driver_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'frame_sha256':{k:v['sha256'] for k,v in x.items()},'exposure':dict(zip(fields,caps)),'mapping':'log1p; nearest voxel display; central two z planes averaged'},indent=2))
