#!/usr/bin/env python3
"""Plot sampled matched-state marginal effects, never interpolated simulation frames."""
import json,pathlib,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=json.loads(pathlib.Path(sys.argv[1]).read_text())
fig,axes=plt.subplots(2,3,figsize=(12,6.4),layout='constrained')
colors=plt.cm.viridis([.08,.25,.42,.59,.76,.93])
for row,seed in enumerate(r['seeds']):
 for col,axis in enumerate((1,2,3)):
  ax=axes[row,col]
  for start,color in zip(r['starts'],colors):
   for sign in (-1,1):
    x=next(x for x in r['results'] if x['seed']==seed and x['start']==start and x['axis']==axis and x['amplitude']==sign*.15)
    ax.plot([t['offset'] for t in x['paired']],[t['delta']['usable'][0] for t in x['paired']],color=color,linestyle='-' if sign==1 else '--',marker='o',ms=2.5,lw=1,label=f't={start:g}' if sign==1 else None)
  
  if max(abs(v) for line in ax.lines for v in line.get_ydata())<1e-7:
   ax.set_ylim(-1e-4,1e-4);ax.text(.5,.85,'No resolved effect at this scale',transform=ax.transAxes,ha='center',fontsize=9)
  ax.axhline(0,color='.65',lw=.6);ax.axvline(2,color='.6',lw=.6,ls=':')
  ax.set_title(f"Founder {seed} · {('investment','allocation','reserve')[col]}")
  ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,2));ax.set_xlabel('Time after decision')
  if col==0:ax.set_ylabel('Producer Δ usable potential')
fig.legend(*axes[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=6,frameon=False,title='Physical checkpoint time · solid +0.15, dashed −0.15 held logit for 0.2')
fig.suptitle('Material actions have thresholded and delayed consequences',fontsize=15)
base=pathlib.Path(sys.argv[2]);fig.savefig(base.with_suffix('.png'),dpi=170);fig.savefig(base.with_suffix('.svg'))
