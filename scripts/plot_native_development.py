#!/usr/bin/env python3
"""Fixed-axis actual-sample figure and compact public readout; no simulation."""
import csv
import hashlib
import html
import io
import json
import pathlib
import subprocess
import sys
from analyze_native_development import read_frame

analysis=pathlib.Path(sys.argv[1]);root=pathlib.Path(sys.argv[2]);out=pathlib.Path(sys.argv[3])
r=json.loads((analysis/'development.json').read_text());index=json.loads((analysis/'manifest-at-analysis.json').read_text())
entries={e['file']:e for e in index['frames']}
frames=[read_frame(root,entries[f['source']['file']]) for f in r['frames']]
svg=['<svg xmlns="http://www.w3.org/2000/svg" width="1500" height="1100" viewBox="0 0 1500 1100">','<rect width="1500" height="1100" fill="white"/>']
def text(x,y,value,size=20,anchor='start'):
    svg.append(f'<text x="{x}" y="{y}" font-family="Arial,Helvetica,sans-serif" font-size="{size}" fill="#263343" text-anchor="{anchor}">{html.escape(str(value))}</text>')
def line(x1,y1,x2,y2,color='#dfe5eb',width=1):
    svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}"/>')
colors=['#1763a0','#db6b25','#3c975e','#a63d81','#7158ab','#66727d','#be9b18','#36a2aa']
text(55,48,'Native garden: sampled development varies; neighbor alignment is mixed',29)
text(55,82,'44 completed actual states, t0–215 · all living IDs · no interpolation or live-state replay',21)
text(55,113,'Cadence 5 = 5,000 physical steps = 6.25 hidden-response times; fast cycles can alias.',21)
csvio=io.StringIO();writer=csv.writer(csvio);writer.writerow(['panel','id','lineage','time','channel','value'])
def panel(left,top,title,ids,field,channels,start,end,labels):
    width=600;height=250
    text(left,top-18,title,22)
    xx=lambda t:left+(t-start)/(end-start)*width
    yy=lambda v:top+(1-v)/2*height
    for v in(-1,0,1):line(left,yy(v),left+width,yy(v));text(left-10,yy(v)+6,v,17,'end')
    for t in ([0,50,100,150,200] if start==0 else [50,65,80,95,110]):
        line(xx(t),top,xx(t),top+height,'#eef1f4');text(xx(t),top+height+25,t,17,'middle')
    for k,(ident,channel,label) in enumerate(zip(ids,channels,labels)):
        actual_field=field[k] if isinstance(field,list) else field
        points=[]
        for f in frames:
            if start<=f['time']<=end+.001 and ident in f['cells']:
                c=f['cells'][ident];v=c[actual_field][channel];points.append((xx(f['time']),yy(v)))
                writer.writerow([title,ident,c['lineages'],f['time'],actual_field+str(channel+1),v])
        svg.append(f'<polyline points="{" ".join(f"{x:.2f},{y:.2f}" for x,y in points)}" fill="none" stroke="{colors[k]}" stroke-width="2"/>')
        for x,y in points:svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.8" fill="{colors[k]}"/>')
        lx=left+(k%4)*155;ly=top+height+58+(k//4)*24
        line(lx,ly-5,lx+15,ly-5,colors[k],3);text(lx+21,ly,label,16)
    text(left+width/2,top+height+112,'Model time · points are observations; lines are guides',17,'middle')
panel(95,190,'Persistent ID 125: hidden channels',[125]*8,'hidden',list(range(8)),0,215,[f'h{i+1}' for i in range(8)])
panel(835,190,'ID 125: adaptation and metabolic actions',[125]*6,['adaptation']*2+['actions']*4,[0,1,3,4,5,6],0,215,['adaptation 1','adaptation 2','import A','import B','growth','respiration'])
panel(95,650,'Selected bonded pair 71 / 154: hidden channel 1',[71,154],'hidden',[0,0],50,110,['ID 71 / lineage 71','ID 154 / lineage 154'])
left=835;top=650;text(left,top-18,'Matched neighbor changes: absolute alignment',22)
for j,w in enumerate(r['neighbor_windows']):
    y=top+j*70;text(left,y+18,w['name'].replace('_',' '),18)
    for k,(name,key,col) in enumerate([('bonded','mean_focal_absolute_bond_cosine','#1763a0'),('nearby unbonded','mean_focal_absolute_control_cosine','#db6b25')]):
        value=w[key];x=left+170;yy=y+k*23
        svg.append(f'<rect x="{x}" y="{yy}" width="{value*600:.2f}" height="16" fill="{col}"/>');text(x+value*600+10,yy+14,f'{value:.3f}',17)
text(left,top+245,'Blue: bonded · Orange: nearby unbonded controls',17)
text(left,top+275,'Same focal ID; matched founder-relatedness and distance.',17)
text(left,top+303,'First differences remove fixed ID/channel offsets.',17)
text(55,1070,'Selected pair alignment is not general evidence of coupling. No autonomous periodicity or evolved organization is established.',20)
svg.append('</svg>');out.parent.mkdir(parents=True,exist_ok=True)
out.with_suffix('.svg').write_text('\n'.join(svg));out.with_suffix('.csv').write_text(csvio.getvalue())
subprocess.run(['rsvg-convert','-o',str(out.with_suffix('.png')),str(out.with_suffix('.svg'))],check=True)
public={k:r[k] for k in ('kind','manifest_sha256','source_commit','pinned_controller_source_sha256','time_bounds','cadence','controller_scales','integrity','limits','frames','neighbor_coverage','neighbor_windows','selected_temporal_pairs','selected_persistent_ids')}
public['windows']=[{k:v for k,v in w.items() if k!='individuals'} for w in r['windows']]
public['selection']='ID125 maximizes observed sign changes among IDs present throughout t0–215; pair71/154 maximizes selected persistent cross-lineage centered first-difference correlation in t50–110. Posthoc examples, not representative estimates.'
public['artifacts']={str(out.with_suffix(s)):hashlib.sha256(out.with_suffix(s).read_bytes()).hexdigest() for s in('.png','.svg','.csv')}
path=out.parent.parent/'native-garden-development.json';path.write_text(json.dumps(public,indent=2))
print(out.with_suffix('.png'))
