#!/usr/bin/env python3
"""Plot recorded native receptor trajectories; no simulation or invented frames."""
import csv, hashlib, html, io, json, pathlib, shutil, subprocess, sys
source=pathlib.Path(sys.argv[1]);out=pathlib.Path(sys.argv[2]);out.parent.mkdir(parents=True,exist_ok=True)
r=json.loads(source.read_text());arms=('off','totals','species');colors={'off':'#687982','totals':'#197c85','species':'#bd5a39'}
svg=['<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="760">','<rect width="1200" height="760" fill="#fbfaf7"/>']
def text(x,y,s,size=17,anchor='start'):
 svg.append(f'<text x="{x}" y="{y}" font-family="Arial,sans-serif" font-size="{size}" fill="#24323b" text-anchor="{anchor}">{html.escape(str(s))}</text>')
def rms(a,b):return (sum((x-y)**2 for x,y in zip(a,b))/len(a))**.5
series={k:[] for k in arms}
for k in arms:
 for a,b in zip(r[k]['trace'],r['off']['trace']):
  assert a['ids']==b['ids'], 'same-ID metric requires matching ordered cohorts'
  series[k].append({'time':a['time'],'hidden_difference_rms':rms(a['hidden'],b['hidden']),
   'action_difference_rms':rms(a['actions'],b['actions']),'bound_total':sum(a['audit']['bound']),'catalyst_total':sum(a['audit']['enzyme'])})
text(35,36,'Giving the controller access to material history',27)
text(35,65,'Native Julia / RX 6750 XT · identical initial cells and fields · three receptor interfaces',18)
text(35,91,'16³ / extent 1³ · dt .02 · 161 observations · the same four cells survive all three windows',16)
data=io.StringIO();writer=csv.writer(data,lineterminator='\n');writer.writerow(['panel','interface','time','value'])
for n,(field,label) in enumerate([('hidden_difference_rms','Hidden-state difference from absent receptors'),('action_difference_rms','Action difference from absent receptors'),('bound_total','Retained food (total physical mass)'),('catalyst_total','Catalyst (total physical mass)')]):
 left=75+(n%2)*590;top=150+(n//2)*280;width=490;height=180
 text(left,top-19,label,18)
 high=max(row[field] for vals in series.values() for row in vals)*1.05 or 1
 xx=lambda t:left+t/32*width;yy=lambda v:top+height-v/high*height
 for f in (0,.5,1):
  y=yy(high*f);svg.append(f'<path d="M {left} {y} h {width}" stroke="#d8dcd9"/>');text(left-8,y+5,f'{high*f:.3g}',13,'end')
 for t in (0,16,32):text(xx(t),top+height+20,t,14,'middle')
 for k in arms:
  points=' '.join(f'{xx(row["time"]):.2f},{yy(row[field]):.2f}' for row in series[k])
  svg.append(f'<polyline points="{points}" fill="none" stroke="{colors[k]}" stroke-width="2.1"/>')
  writer.writerows((field,k,row['time'],row[field]) for row in series[k])
for i,k in enumerate(arms):
 x=75+i*260;svg.append(f'<path d="M {x} 690 h 35" stroke="{colors[k]}" stroke-width="3"/>');text(x+43,695,k,17)
text(35,734,'Different trajectories establish access and coupling; no learning, ecological benefit or converged regime is inferred.',16)
svg.append('</svg>');out.with_suffix('.svg').write_text('\n'.join(svg));out.with_suffix('.csv').write_text(data.getvalue())
if shutil.which('rsvg-convert'):subprocess.run(['rsvg-convert','-o',str(out.with_suffix('.png')),str(out.with_suffix('.svg'))],check=True)
summary={k:r[k] for k in ('gate','timing_off','timing_totals','runtime','source','habitat_source','driver_sha256')}
summary['trace_file_sha256']=hashlib.sha256(source.read_bytes()).hexdigest()
summary['endpoint']={k:{'ids':r[k]['trace'][-1]['ids'],'time':r[k]['trace'][-1]['time'],**series[k][-1], 'residual':r[k]['residual']} for k in arms}
summary['config']={k:r[k]['config'] for k in arms}
summary['timing_overhead_fraction']=r['timing_totals']['median_ms_per_step']/r['timing_off']['median_ms_per_step']-1
out.with_suffix('.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({'endpoint':summary['endpoint'],'timing_overhead_fraction':summary['timing_overhead_fraction']},indent=2))
