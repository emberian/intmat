#!/usr/bin/env python3
"""Verify and compact the source-pinned lattice experiment; no adaptation inference."""
import argparse
import array
import hashlib
import json
from pathlib import Path


def codes(path, count):
    values=array.array('i');values.frombytes(path.read_bytes())
    return {int(x) for x in values[:count] if x>0}


def summarize(directory, benchmark=None):
    directory=Path(directory);runs=json.loads((directory/'results.json').read_text())
    source=hashlib.sha256((directory/'lattice_ecology.jl').read_bytes()).hexdigest()
    rows=[]
    for run in runs:
        assert run['law_source_sha256']==source
        initial,final=run['initial'],run['final'];history=run['history']
        for s in history:
            assert s['finite'] and s['minimum_physical']>=0 and s['hidden_maxabs']<=1.000001
            assert s['capacity_excess']<=2e-6
            assert s['population']==initial['population']+s['births']-s['deaths']
            assert 0<=s['mutations']<=s['births']
        before=next((s for s in history if s['label']=='pre-injury'),None)
        after=next((s for s in history if s['label']=='post-injury'),None)
        if after:assert after['injury_region_population']==0
        path=directory/f"seed-{run['seed']}-{run['arm']}"/'frames'
        if run['arm']=='no_mutation':
            assert final['mutations']==0
            initial_frame=next(path.glob('*-initial.i32'));final_frame=next(path.glob(f"{final['step']:07d}-*.i32"))
            n=1
            for d in run['dimensions']:n*=d
            assert codes(final_frame,n)<=codes(initial_frame,n)
        rows.append(dict(seed=run['seed'],arm=run['arm'],dimensions=run['dimensions'],
            initial_population=initial['population'],population=final['population'],births=final['births'],deaths=final['deaths'],
            hops=final['hops'],mutations=final['mutations'],unique_genotypes=final['unique_genotypes'],
            founder_lineages=final['founder_lineages'],max_generation=final['max_generation'],matrix=final['matrix'],
            injury_population_before=before['injury_region_population'] if before else None,
            injury_population_after=after['injury_region_population'] if after else None,
            injury_population_final=final['injury_region_population'],
            injury_matrix_before=before['injury_region_matrix'] if before else None,
            injury_matrix_final=final['injury_region_matrix'],
            maximum_material_residual=max(abs(s['material_residual']) for s in history),
            maximum_energy_residual=max(abs(s['energy_residual']) for s in history),
            complete_step_seconds=run['step_seconds']/final['step'],voxel_steps_per_second=run['voxel_steps_per_second'],
            resident_array_bytes=final['resident_bytes'],shared_genotype_bytes=final['shared_genotype_bytes'],
            living_topology=final['living_topology'],matrix_topology=final['matrix_topology'],
            final_reserve=final['reserve'],work=final['work']))
    result=dict(model_source_sha256=source,completed_runs=len(rows),all_sampled_states_valid=True,
                maximum_material_residual=max(r['maximum_material_residual'] for r in rows),
                maximum_energy_residual=max(r['maximum_energy_residual'] for r in rows),runs=rows,
                inference='Finite source-pinned mechanism and perturbation study; not a demonstration of adaptive organization.')
    if benchmark:result['complete_step_benchmark']=json.loads(Path(benchmark).read_text())
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory');p.add_argument('--benchmark')
    a=p.parse_args();print(json.dumps(summarize(a.directory,a.benchmark),indent=2))
