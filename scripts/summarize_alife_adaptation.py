#!/usr/bin/env python3
"""Describe preregistered adaptation effects; no significance claims at three seeds."""
import argparse
import json
from pathlib import Path
from statistics import median
import tomllib


def distribution(xs):
    return dict(n=len(xs), minimum=min(xs), median=median(xs), maximum=max(xs)) if xs else dict(n=0)


def summarize(path):
    path = Path(path)
    cfg = tomllib.loads((path / 'protocol.toml').read_text())
    runs = json.loads((path / 'results.json').read_text())
    nulls = json.loads((path / 'life-only.json').read_text())
    for run in runs:
        for assay in run['assays']:
            births, deaths, living = assay['births'], assay['deaths'], assay['living']
            assert [n-b+d for n,b,d in zip(living,births,deaths)] == [cfg['assay_per_group']]*2
    effects = []
    for run in runs:
        effects.append(dict(seed=run['seed'], selection_environment=run['selection_environment'],
                            mutation=run['mutation'], selection_extinct=run['selection_extinct'],
                            selection_population=run['selection_final']['population'],
                            max_living_generation=run['selection_final']['max_generation'],
                            assays=[{k: a[k] for k in ('environment', 'replicate', 'plasticity',
                                    'living', 'births', 'deaths', 'log_frequency_change',
                                    'both_extinct', 'descendant_extinct', 'founder_extinct')}
                                    for a in run['assays']]))
    paired = {}
    for r in runs:
        for a in r['assays']:
            key = (r['seed'], r['selection_environment'], a['environment'], a['replicate'], a['plasticity'])
            paired.setdefault(key, {})[r['mutation']] = a['log_frequency_change']
    deltas = [dict(seed=k[0], selection_environment=k[1], assay_environment=k[2],
                   replicate=k[3], plasticity=k[4], mutation_minus_no_mutation=v[True]-v[False])
              for k, v in paired.items() if True in v and False in v]
    all_assays = [a for r in runs for a in r['assays']]
    all_samples = [json.loads(line) for history in path.glob('**/history.jsonl')
                   for line in history.read_text().splitlines()]
    assert all(s['finite'] for s in all_samples)
    assert all(min(s['source_remaining'][:2]) > 0 for s in all_samples)
    return dict(completed_selection_runs=len(runs), selection_extinctions=sum(r['selection_extinct'] for r in runs),
                completed_assays=len(all_assays), completed_life_only_nulls=len(nulls),
                descendant_effects=distribution([a['log_frequency_change'] for a in all_assays]),
                life_only_effects=distribution([a['log_frequency_change'] for a in nulls]),
                mutation_minus_no_mutation=distribution([d['mutation_minus_no_mutation'] for d in deltas]),
                both_extinct_assays=sum(a['both_extinct'] for a in all_assays),
                both_extinct_life_only_nulls=sum(a['both_extinct'] for a in nulls),
                largest_absolute_material_error=max(abs(s['material_error']) for s in all_samples),
                largest_absolute_chemical_energy_error=max(abs(s['chemical_energy_error']) for s in all_samples),
                effects=effects, paired_deltas=deltas,
                nulls=[{k: a[k] for k in ('seed','environment','replicate','living','births','deaths','log_frequency_change','both_extinct')}
                       for a in nulls])


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_directory')
    args=parser.parse_args()
    print(json.dumps(summarize(args.output_directory), indent=2))
