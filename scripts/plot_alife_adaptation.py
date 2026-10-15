#!/usr/bin/env python3
"""Plot all preregistered competition outcomes and export their exact CSV data."""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COLORS = {'Mutation': '#157eab', 'Founder sorting': '#c47423', 'Same pool': '#737b85'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory')
    parser.add_argument('--output', required=True, help='Output filename stem')
    args = parser.parse_args()
    source, out = Path(args.directory), Path(args.output)
    runs = json.loads((source/'results.json').read_text())
    nulls = json.loads((source/'life-only.json').read_text())
    rows = []
    for r in runs:
        for a in r['assays']:
            rows.append(dict(seed=r['seed'], selection=r['selection_environment'],
                             condition='Mutation' if r['mutation'] else 'Founder sorting',
                             garden=a['environment'], replicate=a['replicate'], plasticity=a['plasticity'],
                             descendants=a['living'][0], ancestors=a['living'][1],
                             descendant_births=a['births'][0], ancestor_births=a['births'][1],
                             descendant_deaths=a['deaths'][0], ancestor_deaths=a['deaths'][1],
                             log_frequency_change=a['log_frequency_change'], both_extinct=a['both_extinct']))
    for a in nulls:
        rows.append(dict(seed=a['seed'], selection='none', condition='Same pool', garden=a['environment'],
                         replicate=a['replicate'], plasticity=a['plasticity'],
                         descendants=a['living'][0], ancestors=a['living'][1],
                         descendant_births=a['births'][0], ancestor_births=a['births'][1],
                         descendant_deaths=a['deaths'][0], ancestor_deaths=a['deaths'][1],
                         log_frequency_change=a['log_frequency_change'], both_extinct=a['both_extinct']))
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.with_suffix('.csv').open('w', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    seeds=sorted({r['seed'] for r in rows})
    fig, axs=plt.subplots(2,len(seeds),figsize=(4*len(seeds),7.2),squeeze=False)
    plt.rcParams.update({'font.family':'DejaVu Sans'})
    count_limit=max(max(r['descendants'],r['ancestors']) for r in rows)*1.08
    conditions=list(COLORS)
    for j,seed in enumerate(seeds):
        share, counts=axs[:,j]
        for x,condition in enumerate(conditions):
            group=[r for r in rows if r['seed']==seed and r['condition']==condition]
            for i,r in enumerate(group):
                # Deterministic visual separation only; no stochastic data modification.
                offset=(i-(len(group)-1)/2)*(.42/max(len(group)-1,1))
                d,a=r['descendants'],r['ancestors']
                marker='o' if r['plasticity']>0 else 's'
                share.scatter(x+offset, d/(d+a) if d+a else -.12,
                              color=COLORS[condition],marker=marker if d+a else 'x',s=34,alpha=.8,
                              linewidths=.8,zorder=3)
                counts.scatter(a,d,color=COLORS[condition],marker=marker,s=34,alpha=.65,linewidths=.6)
            share.text(x,1.075,f'n={len(group)}',ha='center',fontsize=9,color=COLORS[condition])
        failed=[r['selection_environment'] for r in runs if r['seed']==seed and r['selection_extinct']]
        share.set_title(f'Seed {seed}',fontsize=13,fontweight='bold',pad=27)
        share.axhline(.5,color='#999999',lw=.8,ls='--');share.axhline(0,color='#dddddd',lw=.5)
        share.set(xticks=range(3),xticklabels=['Mutation','Founder\nsorting','Same\npool'],ylim=(-.21,1.14),xlim=(-.4,2.4))
        share.set_yticks([-.12,0,.5,1],['Both lost','0','0.5','1'])
        if j==0:share.set_ylabel('Descendant / labeled-group share\nat four lifespans')
        counts.plot([0,count_limit],[0,count_limit],ls='--',color='#aaaaaa',lw=.9)
        counts.set(xlim=(-10,count_limit),ylim=(-10,count_limit),xlabel='Ancestral-genotype group count')
        if j==0:counts.set_ylabel('Descendant-genotype group count')
        counts.set_aspect('equal',adjustable='box')
        if failed:
            share.text(.5,-.31,'Selection extinct: '+', '.join(failed)+' (founder sorting)',
                       transform=share.transAxes,ha='center',fontsize=8,color='#9c5213')
        for ax in (share,counts):
            ax.spines[['right','top']].set_visible(False);ax.tick_params(labelsize=9)
            ax.grid(axis='y',color='#eeeeee',lw=.6,zorder=0)
    fig.suptitle('Inherited descendants compete against their original ancestral genotypes',fontsize=15,y=.985)
    handles=[Line2D([],[],color=c,marker='o',ls='',label=k) for k,c in COLORS.items()]
    handles += [Line2D([],[],color='#444444',marker='o',ls='',label='Plasticity enabled'),
                Line2D([],[],color='#444444',marker='s',ls='',label='Plasticity disabled')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.948),ncol=5,frameon=False,fontsize=9)
    fig.text(.5,.035,'Each point is one assay, not an independent evolutionary seed. Equal starts: 32 + 32; all acquired state reset.\nTwo selection geometries × two gardens × two position repeats; no mutation during assays. Same-pool labels have identical founder pools.\nFailed selection and both-extinct assays are retained. Three evolutionary seeds; descriptive evidence, no significance claim.',
             ha='center',va='center',fontsize=9,linespacing=1.5)
    fig.subplots_adjust(top=.82,bottom=.16,hspace=.46,wspace=.28)
    fig.savefig(out.with_suffix('.png'),dpi=180,facecolor='white')
    fig.savefig(out.with_suffix('.svg'),facecolor='white')
    plt.close(fig)


if __name__=='__main__':
    main()
