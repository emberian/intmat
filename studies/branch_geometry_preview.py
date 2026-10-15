"""Diagnostic projections of exact native graph; no simulated/optical interpolation."""
import json,sys
from PIL import Image, ImageDraw
from pathlib import Path
p=Path(sys.argv[1]);d=json.loads(p.read_text());size=640
im=Image.new('RGB',(size*3,size+60),(8,15,23));draw=ImageDraw.Draw(im)
for panel,(axes,label) in enumerate([((0,1),'XY'),((0,2),'XZ'),((1,2),'YZ')]):
    ax,ay=axes
    def xy(v):return (panel*size+32+(v[ax]-.5)/d['shape'][ax]*(size-64),size-20-(v[ay]-.5)/d['shape'][ay]*(size-64))
    draw.rectangle((panel*size+32,36,(panel+1)*size-32,size-20),outline=(45,60,72))
    for a,b in d['edges']:draw.line((xy(d['positions'][a-1]),xy(d['positions'][b-1])),fill=(156,232,255),width=1)
    draw.text((panel*size+32,12),f"{label}  native {d['shape'][0]}^3  event {d['event']}",fill=(220,230,238))
draw.text((32,size+10),f"Exact face-adjacency: {len(d['positions'])} sites / {len(d['edges'])} edges; {d['tips']} tips / {d['branch_vertices']} branching vertices; bounds {d['bounds']}",fill=(220,230,238))
im.save(p.with_suffix('.png'));print(p.with_suffix('.png'))
