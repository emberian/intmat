"""Plot retained coupled secretion measurements; no interpolated simulation data."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('report',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
r=json.loads(a.report.read_text())['reports']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#f6f3ec','axes.facecolor':'#f6f3ec','axes.edgecolor':'#7e7a73','text.color':'#252b30','axes.labelcolor':'#252b30','xtick.color':'#555c60','ytick.color':'#555c60'})
fig,ax=plt.subplots(1,3,figsize=(14,4.3),gridspec_kw={'width_ratios':[1.15,1,1.1]})
colors=['#157a81','#3553a5','#967025']
for bank,color in zip(r,colors):
 for case in bank['results']:
  d=case['differences'];x=[v['time'] for v in d];y=[v['private'] for v in d]
  ax[0].plot(x,y,color=color,ls='-' if case['raw'] else '--',lw=1.8,alpha=1 if case['raw'] else .65,label=str(bank['seed']) if case['raw'] else None)
  if case['raw']:
   ax[1].plot(x,y,color=color,lw=2)
   t=case['private_payback']['first'];value=next(v['private'] for v in d if v['time']==t);ax[1].scatter([t],[value],color=color,s=35,zorder=3)
   ax[1].annotate(f'{t:g}',(t,value),xytext=(0,9),textcoords='offset points',ha='center',color=color)
case=r[0]['results'][1];d=case['pulse'];x=[v['time'] for v in d]
for key,label,color in [('solubilized','Raw food released','#157a81'),('enzyme','Catalyst remaining','#ba654b'),('bound','Food retained in matrix','#77569c')]:
 ax[2].plot(x,[sum(v[key]) for v in d],label=label,color=color,lw=2)
for axis in ax[:2]:
 axis.axhline(0,color='#777d80',lw=.8);axis.set_ylabel('Producer usable potential\nchange versus no secretion')
 axis.axvspan(0,.5,color='#c2b597',alpha=.2,lw=0)
for axis in ax:axis.set_xlabel('Physical time');axis.grid(axis='y',alpha=.15)
ax[0].set_title('Construction pays when raw food is present',loc='left',fontsize=11,pad=12)
ax[0].set_xlim(0,16);ax[0].legend(title='Founder bank',frameon=True,facecolor='#f6f3ec',edgecolor='none',fontsize=9,loc='lower right')
ax[0].text(.03,.97,'Solid: raw food available\nDashed: no raw food',transform=ax[0].transAxes,va='top',fontsize=9)
ax[1].set_title('First positive samples: 0.4–0.6',loc='left',fontsize=11,pad=12);ax[1].set_xlim(0,1.2)
ax[2].set_title('The pulse leaves material behind',loc='left',fontsize=11,pad=12);ax[2].set_ylabel('Material amount · founder 211');ax[2].legend(frameon=False,fontsize=8,loc='center right')
fig.suptitle('A paid secretion pulse opens a local food pocket',x=.06,ha='left',fontsize=17,y=1.02)
fig.text(.06,-.03,'Programmed investment 0.8 for 0.5 time (shaded); production limits actual secretion by available stores.\nThree matched four-cell experiments, observed through 16 time units. This is a physical opportunity assay, not learned cooperation.',fontsize=9,color='#555c60')
fig.tight_layout();a.output.parent.mkdir(parents=True,exist_ok=True)
for ext in ('svg','png'):fig.savefig(a.output.with_suffix('.'+ext),dpi=160,bbox_inches='tight')
