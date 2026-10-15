#!/usr/bin/env python3
"""Plot sampled recorded metrics and original-source-restored archived positions."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--summary',type=Path,default=Path('docs/cellular-soak-completion.json'))
p.add_argument('--output',type=Path,default=Path('docs/figures'))
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
s=json.loads(a.summary.read_text());frames=s['archived_visual_frames']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
for run,color in zip(s['runs'],['#246b9a','#b85c26']):
 rows=run['sampled_trajectory'];name=Path(run['run']).name
 x=[(r['observed_at_unix']-rows[0]['observed_at_unix'])/3600 for r in rows]
 for ax,key,title in zip(axes.flat,['population','bonds','surviving_lineages','max_generation'],
                         ['Living cells','Recorded bonds','Surviving founder lineages','Maximum living generation']):
  ax.plot(x,[r['metrics'][key] for r in rows],color=color,label=name,lw=1.8)
  ax.set(xlabel='Elapsed wall time (hours)',title=title,xlim=(0,12));ax.grid(alpha=.18)
axes[0,0].legend(frameon=False);axes[1,0].set_yscale('log');axes[1,0].set_ylim(.8,80)
fig.suptitle('Two original cellular soaks completed their 12-hour wall deadline',fontsize=15)
fig.supxlabel('Recorded metrics sampled every ≈10 wall minutes; narrow events can fall between plotted points',fontsize=10)
fig.savefig(a.output/'cellular-soak-completion-metrics.png',dpi=150);plt.close(fig)
fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
for ax,frame in zip(axes.flat,frames):
 ex,ey,_=frame['grid']['extent'];generation_max=max([agent['generation'] for agent in frame['recorded_positions']] or [1])
 for agent in frame['recorded_positions']:
  x,y,_=agent['position'];radius=agent['radius']
  ax.add_patch(Circle((x,y),radius,facecolor=plt.cm.viridis(agent['generation']/max(1,generation_max)),edgecolor='white',linewidth=.3,alpha=.9))
 ax.set(xlim=(-ex/2,ex/2),ylim=(-ey/2,ey/2),aspect='equal',
        title=f"{frame['run']} · {frame['label']}\nt={frame['time']:,.0f}; N={len(frame['recorded_positions'])}; bonds={frame['bonds']}")
 ax.set_xticks([]);ax.set_yticks([])
fig.suptitle('Actual archived cell positions, restored with the original source',fontsize=15)
fig.supxlabel('Circle radius is recorded cell radius. Color shows generation within each panel; no bonds are inferred from proximity.',fontsize=10)
fig.savefig(a.output/'cellular-soak-completion-positions.png',dpi=150);plt.close(fig)
