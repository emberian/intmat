#!/usr/bin/env python3
"""Draw measured catalytic fields, physical disks and persistent bonds.

No interpolation in time. Exposure is fixed across the supplied recorded frames.
Requires Pillow. Field interpolation is nearest-neighbour to expose grid support.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('source',type=Path);p.add_argument('output',type=Path)
p.add_argument('--panel',type=int,default=440)
a=p.parse_args();paths=sorted(a.source.glob('frame-*.json'));assert paths
states=[json.loads(path.read_text()) for path in paths]
a.output.mkdir(parents=True,exist_ok=True)
keys=[('enzyme_a','enzyme_b'),('precursor_a','precursor_b'),('retained_a','retained_b')]
titles=['Constructed catalyst', 'Accessible food', 'Retained food', 'Matrix transport resistance']
exposure=[]
for pair in keys:
    exposure.append(max(max(x+y for x,y in zip(s['fields'][pair[0]]['values'],s['fields'][pair[1]]['values'])) for s in states))
font_path=next((p for p in ('/System/Library/Fonts/Supplemental/Arial.ttf','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf') if Path(p).exists()),None)
def make_font(size):
    return ImageFont.truetype(font_path,size) if font_path else ImageFont.load_default(size=size)
font=make_font(16);small=make_font(12);big=make_font(22)
for number,(path,s) in enumerate(zip(paths,states)):
    width=a.panel;canvas=Image.new('RGB',(2*width+48,2*width+176),(13,20,24));draw=ImageDraw.Draw(canvas)
    draw.text((16,12),'A place that holds its chemistry',font=big,fill=(239,231,213))
    m=s['metrics'];cat=m['catalytic']
    draw.text((16,43),f"t {s['time']:.1f}   ·   {m['population']} living cells   ·   {m['bonds']} bonds   ·   {m['births']} births / {m['deaths']} deaths",font=font,fill=(190,201,196))
    lx,ly,_=s['grid']['extent'];centers={c['id']:c['position'][:2] for c in s['agents']}
    for panel in range(4):
        pair=keys[min(panel,2)]
        nx,ny=s['fields'][pair[0]]['shape'];pixels=[]
        for index,(av,bv) in enumerate(zip(s['fields'][pair[0]]['values'],s['fields'][pair[1]]['values'])):
            total=av+bv;mix=av/total if total else .5
            brightness=math.log1p(20*total/max(exposure[min(panel,2)],1e-30))/math.log(21)
            amber=(235,153,49);teal=(40,193,185)
            if panel==3:
                resistance=1-s['fields']['matrix_permeability']['values'][index]
                pixels.append(tuple(round(12+resistance*v) for v in (196,169,221)))
            else:
                pixels.append(tuple(round(12+(mix*amber[k]+(1-mix)*teal[k])*brightness*.9) for k in range(3)))
        image=Image.new('RGB',(nx,ny));image.putdata(pixels)
        image=image.transpose(Image.Transpose.FLIP_TOP_BOTTOM).resize((width,width),Image.Resampling.NEAREST)
        layer=Image.new('RGBA',image.size);d=ImageDraw.Draw(layer)
        def point(x,y):return ((x/lx+.5)*width,(.5-y/ly)*width)
        for ia,ib in m['cell_bonds']:
            if ia not in centers or ib not in centers:continue
            x,y=centers[ia];xx,yy=centers[ib]
            dx=(xx-x+lx/2)%lx-lx/2;dy=(yy-y+ly/2)%ly-ly/2
            for ox in (-lx,0,lx):
                for oy in (-ly,0,ly):d.line([point(x+ox,y+oy),point(x+ox+dx,y+oy+dy)],fill=(235,245,224,140),width=1)
        for c in s['agents']:
            x,y=c['position'][:2];rx=c['radius']/lx*width;ry=c['radius']/ly*width
            for ox in (-lx,0,lx):
                for oy in (-ly,0,ly):
                    px,py=point(x+ox,y+oy)
                    if px+rx<0 or py+ry<0 or px-rx>width or py-ry>width:continue
                    d.ellipse((px-rx,py-ry,px+rx,py+ry),fill=(7,12,15,90),outline=(245,241,220,220),width=1)
        image=Image.alpha_composite(image.convert('RGBA'),layer).convert('RGB')
        x=16+(panel%2)*(width+16);y=87+(panel//2)*(width+28);canvas.paste(image,(x,y))
        draw.text((x,y-20),titles[panel],font=font,fill=(235,228,213))
    draw.text((16,2*width+122),f"Enzyme mass {sum(cat['enzyme_mass']):.3f}   ·   Retained {sum(m['matrix']['bound_mass']):.3f}   ·   A amber / B teal",font=font,fill=(205,214,204))
    draw.text((16,2*width+148),'Recorded fields and physical disks. Fixed log exposures; resistance uses fixed 0–1 scale.',font=small,fill=(140,159,155))
    canvas.save(a.output/f'frame-{number:04d}.png')
metadata=dict(source=str(a.source.resolve()),frames=len(paths),first_time=states[0]['time'],last_time=states[-1]['time'],
    exposure=exposure,source_sha256=[hashlib.sha256(path.read_bytes()).hexdigest() for path in paths],
    mapping='Exact sampled field values; fixed log1p exposure; nearest-neighbour field display. Physical radii/positions and recorded periodic bonds; no time interpolation.')
(a.output/'manifest.json').write_text(json.dumps(metadata,indent=2))
print(json.dumps(dict(frames=len(paths),output=str(a.output))))
