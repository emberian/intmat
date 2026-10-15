#!/usr/bin/env python3
"""Derived neural/neighbor readout from completed native observatory exports only."""
import argparse
import array
import collections
import hashlib
import json
import math
import pathlib
import statistics
import subprocess
import sys


def mean(values):
    return statistics.fmean(values) if values else None


def rms_std(rows, field):
    return math.sqrt(statistics.fmean(statistics.variance(r[field][j] for r in rows)
        for j in range(len(rows[0][field])))) if len(rows) > 1 else 0.


def dot_cos(a, b):
    aa = sum(v*v for v in a);bb = sum(v*v for v in b)
    return sum(x*y for x, y in zip(a, b))/math.sqrt(aa*bb) if aa*bb > 1e-24 else None


def pooled_corr(a, b):
    if len(a) < 3:
        return None
    ma = [statistics.fmean(row[j] for row in a) for j in range(len(a[0]))]
    mb = [statistics.fmean(row[j] for row in b) for j in range(len(b[0]))]
    return dot_cos([v-ma[j] for row in a for j,v in enumerate(row)],
                   [v-mb[j] for row in b for j,v in enumerate(row)])


def read_frame(root, entry):
    path = root/entry['file'];raw = path.read_bytes();meta = json.loads(raw)
    binary = (path.parent/meta['binary']).read_bytes()
    assert hashlib.sha256(binary).hexdigest() == meta['binary_sha256']
    assert meta['source_sha256'] == entry['source_sha256']
    def values(name):
        spec=meta['arrays'][name];a=array.array({'i32':'i','f32':'f','f64':'d'}[spec['type']])
        a.frombytes(binary[spec['offset']:spec['offset']+spec['bytes']])
        if meta['endian'] != sys.byteorder:a.byteswap()
        return a
    ids=list(values('ids'));n=len(ids);decoded={name:values(name) for name in
        ('parents','lineages','hidden','adaptation','actions','energy','position','radii','mass','metabolites','last_step_secretion')}
    cells={ident: {'id':ident,'time':meta['time'],'step':meta['step'],
        **{name:[decoded[name][i+n*j] for j in range(meta['arrays'][name]['shape'][1])]
           for name in ('hidden','adaptation','actions','position','metabolites','last_step_secretion')},
        **{name:decoded[name][i] for name in ('parents','lineages','energy','radii','mass')}} for i,ident in enumerate(ids)}
    slots=values('bond_ids');directed={(ids[i],slots[i+n*k]) for i in range(n)
        for k in range(meta['arrays']['bond_ids']['shape'][1]) if slots[i+n*k]>0}
    edges={tuple(sorted((u,v))) for u,v in directed if u!=v and v in cells and (v,u) in directed}
    for c in cells.values():
        exp=[math.exp(2*v) for v in c['actions'][3:7]];c['metabolic']=[v/sum(exp) for v in exp]
    return {'time':meta['time'],'step':meta['step'],'extent':meta['extent'],'cells':cells,'edges':edges,
        'source':{'file':entry['file'],'metadata_sha256':hashlib.sha256(raw).hexdigest(),
            'binary_sha256':meta['binary_sha256'],'source_sha256':meta['source_sha256']}}


def dist_ratio(a,b,extent):
    ds=[abs(a['position'][j]-b['position'][j]) for j in range(3)]
    return math.sqrt(sum(min(v,L-v)**2 for v,L in zip(ds,extent)))/(a['radii']+b['radii'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('observatory',type=pathlib.Path)
    p.add_argument('output',type=pathlib.Path);p.add_argument('--max-time',type=float,default=215.001)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    raw=(args.observatory/'observatory.json').read_bytes();index=json.loads(raw)
    (args.output/'manifest-at-analysis.json').write_bytes(raw)
    grouped=collections.defaultdict(list)
    for entry in index['frames']:
        if entry['time']<=args.max_time:grouped[entry['step']].append(entry)
    frames=[];duplicates=0
    for step,entries in sorted(grouped.items()):
        chosen=max(entries,key=lambda e:e['side']);frame=read_frame(args.observatory,chosen)
        for other in entries:
            if other is not chosen:
                duplicate=read_frame(args.observatory,other)
                assert duplicate['cells']==frame['cells'] and duplicate['edges']==frame['edges']
                duplicates+=1
        frames.append(frame)
    cadence=statistics.median(b['time']-a['time'] for a,b in zip(frames,frames[1:]))
    assert all(abs(b['time']-a['time']-cadence)<1e-5 for a,b in zip(frames,frames[1:]))
    by_id=collections.defaultdict(list)
    for f in frames:
        for ident,c in f['cells'].items():by_id[ident].append(c)
    window_bounds=[('startup',0.,50.001),('pre_decline',50.,110.001),('decline',110.,frames[-1]['time']+.001)]
    windows=[]
    for name,start,end in window_bounds:
        selected=[f for f in frames if start<=f['time']<=end];individuals=[]
        for ident,rows in by_id.items():
            rows=[r for r in rows if start<=r['time']<=end]
            if len(rows)<6:continue
            pairs=[(a,b) for a,b in zip(rows,rows[1:]) if b['step']-a['step']==5000]
            individuals.append({'id':ident,'parent':rows[0]['parents'],'lineage':rows[0]['lineages'],
                'samples':len(rows),'time_bounds':[rows[0]['time'],rows[-1]['time']],
                'complete_window':len(rows)==len(selected),
                'rms_channel_std':{field:rms_std(rows,field) for field in ('hidden','adaptation','actions','metabolic')},
                'hidden_saturated_fraction':sum(abs(v)>.95 for r in rows for v in r['hidden'])/(8*len(rows)),
                'hidden_sign_changes':sum((x<0)!=(y<0) for a,b in pairs for x,y in zip(a['hidden'],b['hidden'])),
                'mean_rms_hidden_change':mean([math.sqrt(sum((x-y)**2 for x,y in zip(a['hidden'],b['hidden']))/8) for a,b in pairs])})
        windows.append({'name':name,'bounds':[start,end],'frames':len(selected),'qualifying_ids':len(individuals),
            'complete_window_ids':sum(x['complete_window'] for x in individuals),
            'median_within_id_std':{field:statistics.median(x['rms_channel_std'][field] for x in individuals) if individuals else None for field in ('hidden','adaptation','actions','metabolic')},
            'median_saturated_fraction':statistics.median(x['hidden_saturated_fraction'] for x in individuals) if individuals else None,
            'ids_with_sign_changes':sum(x['hidden_sign_changes']>0 for x in individuals),'individuals':individuals})
    comparisons=[];coverage=collections.Counter();coverage_intervals=[]
    for old,new in zip(frames,frames[1:]):
        common=old['cells'].keys()&new['cells'].keys()
        delta={i:[b-a for a,b in zip(old['cells'][i]['hidden'],new['cells'][i]['hidden'])] for i in common}
        active={i for i,v in delta.items() if math.sqrt(sum(x*x for x in v)/8)>=.01}
        endpoint_edges=old['edges']&new['edges'];excluded=old['edges']|new['edges']
        coverage_intervals.append({'midtime':(old['time']+new['time'])/2,
            'active_directed_endpoint_bonds':sum(u in active and v in active for u,v in endpoint_edges)*2})
        for u,v in endpoint_edges:
            for focal,partner in ((u,v),(v,u)):
                coverage['directed_endpoint_bonds']+=1
                if focal not in active or partner not in active:continue
                coverage['active_directed_endpoint_bonds']+=1
                same=new['cells'][focal]['lineages']==new['cells'][partner]['lineages']
                target=dist_ratio(new['cells'][focal],new['cells'][partner],new['extent'])
                candidates=[]
                for control in active-{focal,partner}:
                    if tuple(sorted((focal,control))) in excluded:continue
                    if (new['cells'][focal]['lineages']==new['cells'][control]['lineages'])!=same:continue
                    ratio=dist_ratio(new['cells'][focal],new['cells'][control],new['extent'])
                    if abs(ratio-target)<=.25:candidates.append((abs(ratio-target),control,ratio))
                if not candidates:continue
                gap,control,ratio=min(candidates)
                comparisons.append({'time':new['time'],'focal':focal,'partner':partner,'control':control,
                    'same_founder_lineage':same,'bond_distance_radii':target,'control_distance_radii':ratio,'distance_gap':gap,
                    'bond_delta_cosine':dot_cos(delta[focal],delta[partner]),
                    'control_delta_cosine':dot_cos(delta[focal],delta[control])})
    neighbor_windows=[]
    for name,start,end in window_bounds:
        matches=[x for x in comparisons if start<=x['time']-cadence/2<end]
        focal=collections.defaultdict(list)
        for x in matches:focal[(x['time'],x['focal'])].append(x)
        neighbor_windows.append({'name':name,'directed_matched_comparisons':len(matches),'focal_time_observations':len(focal),
            'active_directed_endpoint_bonds':sum(x['active_directed_endpoint_bonds'] for x in coverage_intervals if start<=x['midtime']<end),
            'same_lineage_comparisons':sum(x['same_founder_lineage'] for x in matches),
            'mean_focal_bond_delta_cosine':mean([statistics.fmean(x['bond_delta_cosine'] for x in row) for row in focal.values()]),
            'mean_focal_control_delta_cosine':mean([statistics.fmean(x['control_delta_cosine'] for x in row) for row in focal.values()]),
            'mean_focal_excess':mean([statistics.fmean(x['bond_delta_cosine']-x['control_delta_cosine'] for x in row) for row in focal.values()]),
            'mean_focal_absolute_bond_cosine':mean([statistics.fmean(abs(x['bond_delta_cosine']) for x in row) for row in focal.values()]),
            'mean_focal_absolute_control_cosine':mean([statistics.fmean(abs(x['control_delta_cosine']) for x in row) for row in focal.values()]),
            'mean_focal_absolute_excess':mean([statistics.fmean(abs(x['bond_delta_cosine'])-abs(x['control_delta_cosine']) for x in row) for row in focal.values()]),
            'median_distance_match_gap':statistics.median(x['distance_gap'] for x in matches) if matches else None})
    # Selected persistent-pair temporal diagnostic, not a permutation p-value.
    temporal=[]
    for name,start,end in window_bounds[1:]:
        fs=[f for f in frames if start<=f['time']<=end];present=set.intersection(*(set(f['cells']) for f in fs))
        counts=collections.Counter(edge for f in fs for edge in f['edges'])
        for (u,v),count in counts.items():
            if count<6 or u not in present or v not in present:continue
            if fs[0]['cells'][u]['lineages']==fs[0]['cells'][v]['lineages']:continue
            a=[[y-x for x,y in zip(f0['cells'][u]['hidden'],f1['cells'][u]['hidden'])] for f0,f1 in zip(fs,fs[1:])]
            b=[[y-x for x,y in zip(f0['cells'][v]['hidden'],f1['cells'][v]['hidden'])] for f0,f1 in zip(fs,fs[1:])]
            observed=pooled_corr(a,b)
            if observed is None:continue
            shifts=[{'offset':k,'time_offset':k*cadence,'correlation':pooled_corr(a,b[k:]+b[:k])} for k in range(1,len(b))]
            temporal.append({'window':name,'ids':[u,v],'bonded_frames':count,'window_frames':len(fs),
                'lineages':[fs[0]['cells'][u]['lineages'],fs[0]['cells'][v]['lineages']],
                'centered_difference_correlation':observed,'circular_shift_diagnostic':shifts})
    temporal=[p for name,_,_ in window_bounds[1:] for p in sorted(
        (x for x in temporal if x['window']==name),key=lambda x:abs(x['centered_difference_correlation']),reverse=True)[:4]]
    persistent=[i for i,rows in by_id.items() if len(rows)==len(frames)]
    ranking=[]
    for i in persistent:
        rows=[r for r in by_id[i] if r['time']>=30]
        ranking.append((sum((a<0)!=(b<0) for x,y in zip(rows,rows[1:]) for a,b in zip(x['hidden'],y['hidden'])),rms_std(rows,'hidden'),i))
    ranked=sorted(ranking,reverse=True);chosen_ids=[ranked[0][2],min(ranking)[2]] if ranking else []
    frame_summary=[{'time':f['time'],'step':f['step'],'population':len(f['cells']),'edges':len(f['edges']),
        'mean_energy':mean([c['energy'] for c in f['cells'].values()]),
        'hidden_saturated_fraction':sum(abs(v)>.95 for c in f['cells'].values() for v in c['hidden'])/(8*len(f['cells'])),
        'source':f['source']} for f in frames]
    law=subprocess.check_output(['git','show',index['source_commit']+':src/device_cellular_kernels.jl'])
    result={'kind':'read-only derived native neural development and neighbor relationships',
        'manifest_sha256':hashlib.sha256(raw).hexdigest(),'source_commit':index['source_commit'],
        'pinned_controller_source_sha256':hashlib.sha256(law).hexdigest(),
        'time_bounds':[frames[0]['time'],frames[-1]['time']],'cadence':cadence,
        'controller_scales':{'physical_dt':index['physical_dt'],'response_time':.8,'adaptation_time':8.,'adaptation_strength':1.5,'plasticity_rate':.002,'plasticity_decay':.02,'cadence_over_response_time':cadence/.8,'cadence_over_adaptation_time':cadence/8},
        'integrity':{'selected_frames':len(frames),'equal_cell_state_duplicate_exports':duplicates,'all_metadata_and_binary_hashes_verified':True},
        'limits':['time5 cadence is5000physicalsteps and6.25response times; periods below10 can alias and apparent alternation is not proof of periodicity',
            'within-ID differences remove fixed offsets but not inherited dynamics, common chemical forcing, division or selection',
            'neighbor controls share focal cell and founder-relatedness class, with normalized-distance gap<=.25; conditioning on RMS hidden change>=.01',
            'endpoint bonds need not persist continuously between observations; comparisons are dependent descriptive units, not independent trials',
            'circular shifts preserve each ID difference sequence; this is an alignment diagnostic, not a significance test or proof of causal coupling',
            'all observed cells retained; statistics requiring>=6samples and complete cohorts are explicitly survivor-conditioned; no adaptation/evolved organization claim'],
        'frames':frame_summary,'windows':windows,'neighbor_coverage':dict(coverage),'neighbor_windows':neighbor_windows,
        'selected_temporal_pairs':temporal,'selected_persistent_ids':chosen_ids,
        'selected_actual_traces':{str(i):by_id[i] for i in chosen_ids}}
    (args.output/'development.json').write_text(json.dumps(result,indent=2))
    (args.output/'neighbor-matches.json').write_text(json.dumps(comparisons,separators=(',',':')))
    print(json.dumps({'output':str(args.output),'frames':len(frames),'matched_neighbor_comparisons':len(comparisons),'selected_ids':chosen_ids}))


if __name__=='__main__':main()
