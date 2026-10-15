#!/usr/bin/env python3
"""Render actual Julia assay outputs with standard-library SVG; no simulation."""
import csv, hashlib, html, io, json, pathlib, shutil, subprocess, sys
root=pathlib.Path(sys.argv[1]);out=pathlib.Path(sys.argv[2]);out.parent.mkdir(parents=True,exist_ok=True)
r=json.loads((root/'report.json').read_text())
arms={k:json.loads((root/(k+'.json')).read_text()) for k in r['final']}
svg=['<svg xmlns="http://www.w3.org/2000/svg" width="1260" height="900" viewBox="0 0 1260 900">','<rect width="1260" height="900" fill="#fbfaf7"/>']
def text(x,y,s,size=18,anchor='start'):
    svg.append(f'<text x="{x}" y="{y}" font-family="Arial,sans-serif" font-size="{size}" fill="#24323b" text-anchor="{anchor}">{html.escape(str(s))}</text>')
colors={'none':'#637783','neural':'#9f3455','marks':'#207c91','retention':'#bb722c'}
csvout=io.StringIO();writer=csv.writer(csvout, lineterminator="\n");writer.writerow(['panel','condition','time','value'])
text(40,40,'Native garden: memory changes response to accessible food',27)
text(40,72,f"Common state at t{r['warm']:g}; four founders; native 16³ / 1³ world; dt {r['dt']:.3g}, observation .2",18)
text(40,100,'Erasure × pulse interaction is a coupled causal response, not a score for learning or adaptation.',17)
def panel(left,top,title,series):
    width,height=500,245;text(left,top-20,title,21)
    allv=[v for pairs in series.values() for _,v in pairs if v is not None]
    low=min(0,min(allv));high=max(allv);high=high if high>low else low+1
    pad=.04*(high-low);low-=pad;high+=pad
    mint=r['warm'];maxt=mint+r['horizon'];xx=lambda t:left+(t-mint)/(maxt-mint)*width
    yy=lambda v:top+height-(v-low)/(high-low)*height
    for fraction in (0,.5,1):
        v=low+(high-low)*fraction;y=yy(v)
        svg.append(f'<path d="M {left} {y} h {width}" stroke="#d8dcd9"/>');text(left-8,y+5,f'{v:.3g}',14,'end')
    for t in (mint,mint+r['horizon']/2,maxt):text(xx(t),top+height+23,f'{t:g}',15,'middle')
    for n,(name,pairs) in enumerate(series.items()):
        points=' '.join(f'{xx(t):.2f},{yy(v):.2f}' for t,v in pairs if v is not None)
        svg.append(f'<polyline points="{points}" fill="none" stroke="{colors[name]}" stroke-width="2.2"/>')
        text(left+n*125,top+height+54,name,16)
        svg.append(f'<path d="M {left+n*125} {top+height+63} h 65" stroke="{colors[name]}" stroke-width="3"/>')
        writer.writerows((title,name,t,v) for t,v in pairs)
for col,(field,title) in enumerate([('hidden_interaction_rms','Hidden-state pulse interaction (RMS)'),('action_interaction_rms','Action pulse interaction (RMS)')]):
    panel(90+col*610,160,title,{k:[(x['time'],x[field]) for x in rows] for k,rows in r['interactions'].items()})
for col,(field,title) in enumerate([('bound','Bound food: pulse − sham (mass)'),('catalyst','Catalyst: pulse − sham (mass)')]):
    series={}
    for name in colors:
        pairs=zip(arms[name+'_pulse']['trace'],arms[name+'_sham']['trace'])
        series[name]=[(a['time'],a[field]-b[field]) for a,b in pairs]
    panel(90+col*610,550,title,series)
text(40,880,'Same-ID intersection is retained at every sample. Lines connect observations. This small periodic world is a bounded mechanism assay.',16)
svg.append('</svg>');out.with_suffix('.svg').write_text('\n'.join(svg));out.with_suffix('.csv').write_text(csvout.getvalue())
if shutil.which('rsvg-convert'):subprocess.run(['rsvg-convert','-o',str(out.with_suffix('.png')),str(out.with_suffix('.svg'))],check=True)
summary={k:v for k,v in r.items() if k!='interactions'}
summary['interactions']={k:{'peak_hidden_rms':max(x['hidden_interaction_rms'] or 0 for x in v),'final_hidden_rms':v[-1]['hidden_interaction_rms'],'peak_action_rms':max(x['action_interaction_rms'] or 0 for x in v),'final_action_rms':v[-1]['action_interaction_rms'],'minimum_common_ids':min(len(x['common_ids']) for x in v)} for k,v in r['interactions'].items()}
for name, rows in r['interactions'].items():
    for field in ('hidden_interaction_rms','action_interaction_rms'):
        peak=max(rows,key=lambda x:x[field] or 0)
        summary['interactions'][name][field+'_peak_delay']=peak['time']-r['warm']
        summary['interactions'][name][field+'_rms_over_window']=(sum((x[field] or 0)**2 for x in rows)/len(rows))**.5
summary['files']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in root.glob('*.json')}
summary['branch_endpoint']={k:{name:v['trace'][-1][name] for name in ('population','births','deaths','free','raw','bound','catalyst','permeability_mean','mark_rms')} for k,v in arms.items()}
out.with_suffix('.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary['interactions'],indent=2))
