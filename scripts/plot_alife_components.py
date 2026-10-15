#!/usr/bin/env python3
"""Plot exploratory ancestry-matched component outcomes, preserving every assay."""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

COMPONENTS=['evolved_both','ancestral_controller','ancestral_metabolism','ancestral_both']
LABELS=['Evolved controller + metabolism','Ancestral controller; evolved metabolism',
        'Evolved controller; ancestral metabolism','Ancestral controller + metabolism']


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory');p.add_argument('--output',required=True)
    args=p.parse_args();source=Path(args.directory);out=Path(args.output)
    results=json.loads((source/'results.json').read_text());seeds=sorted({r['seed'] for r in results})
    values=np.full((8,len(seeds)),np.nan);lookup={}
    rows=[]
    for r in results:
        d,a=r['living'];row=2*COMPONENTS.index(r['component'])+r['replicate']-1;col=seeds.index(r['seed'])
        if d+a:values[row,col]=d/(d+a)
        lookup[row,col]=r
        rows.append({k:r[k] for k in ['seed','replicate','component','log_frequency_change','both_extinct']} |
                    dict(test=d,ancestor=a,test_births=r['births'][0],ancestor_births=r['births'][1],
                         test_deaths=r['deaths'][0],ancestor_deaths=r['deaths'][1]))
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.with_suffix('.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    fig,ax=plt.subplots(figsize=(10,6.7));cmap=LinearSegmentedColormap.from_list('shares',['#c47423','#f2f1ed','#157eab'])
    cmap.set_bad('#dddddd');im=ax.imshow(values,vmin=0,vmax=1,cmap=cmap,aspect='auto')
    ax.set_xticks(range(len(seeds)),[f'Seed {s}' for s in seeds]);ax.xaxis.tick_top();ax.tick_params(axis='both',length=0,pad=9)
    ax.set_yticks(range(8),[f'{label} · repeat {rep}' for label in LABELS for rep in (1,2)],fontsize=10)
    for row in range(8):
        for col in range(len(seeds)):
            r=lookup.get((row,col));
            if not r:continue
            d,a=r['living'];share=values[row,col]
            ax.text(col,row,f'{d} / {a}',ha='center',va='center',fontsize=13,
                    color='white' if share>.8 or share<.15 else '#222222')
    for y in (1.5,3.5,5.5):ax.axhline(y,color='white',lw=3)
    for spine in ax.spines.values():spine.set_visible(False)
    fig.suptitle('Which inherited component retains performance?',fontsize=15,y=.98)
    fig.text(.68,.9,'Terminal test / ancestry-matched ancestor counts',ha='center',fontsize=10)
    bar=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.055,pad=.1,ticks=[0,.5,1])
    bar.set_label('Test-group share; gray cells mean both groups extinct',fontsize=10)
    fig.text(.5,.035,'Exploratory, separate from the fixed primary study. Each cell is one assay; three evolutionary seeds.\nEqual 32 + 32 starts; uniform garden; plasticity and mutation disabled; acquired state reset.\nExact ancestral genotypes, sampling and positions matched across components. No significance claim.',
             ha='center',va='center',fontsize=9,linespacing=1.4)
    fig.subplots_adjust(left=.47,right=.97,top=.85,bottom=.2)
    fig.savefig(out.with_suffix('.png'),dpi=180,facecolor='white');fig.savefig(out.with_suffix('.svg'),facecolor='white')
    plt.close(fig)


if __name__=='__main__':main()
