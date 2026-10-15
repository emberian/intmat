"""Check saved repair histories independently of the Julia producer.

This verifies persistent identities and accounting across actual captured states;
it does not infer unrecorded motion or certify continuum convergence.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def check(run_directory):
    root = Path(run_directory)
    run = json.loads((root / 'run.json').read_text())
    previous = {}
    checks = 0
    maximum_stock_error = 0.0
    maximum_energy_fraction = 0.0
    sources = []
    for frame, history in zip(run['frames'], run['history'], strict=True):
        path = root / frame['directory'] / 'geometry.json'
        g = json.loads(path.read_text())
        n = len(g['bond_ids'])
        assert g['array_order'] == 'F'
        assert len(set(g['bond_ids'])) == n
        assert g['bond_ids'] == sorted(g['bond_ids'])
        assert len(g['positions']) == 3 * len(g['node_ids'])
        assert all(math.isfinite(x) for x in g['positions'] + g['velocities'])
        assert min(g['repair_stock']) >= 0
        assert len(g['bonds']) == 2 * n
        assert len(g['parent_bond']) == len(g['generation']) == n
        current = {}
        for i, identity in enumerate(g['bond_ids']):
            edge = (g['bonds'][i], g['bonds'][n+i])
            damage, parent, generation = g['damage'][i], g['parent_bond'][i], g['generation'][i]
            assert damage in (0, 1)
            assert g['rest_length'][i] > 0
            if identity in previous:
                old = previous[identity]
                assert edge == old['edge']
                assert damage >= old['damage']
                assert g['rest_length'][i] == old['rest']
                assert parent == old['parent']
            if parent:
                ancestor = current[parent]
                assert ancestor['damage'] == 1
                assert edge == ancestor['edge']
                assert generation == ancestor['generation'] + 1 <= 2
                assert g['birth_time'][i] >= ancestor['birth']
                assert math.isclose(g['stiffness'][i], .7 * ancestor['stiffness'], rel_tol=1e-14)
            else:
                assert generation == 0 and g['birth_time'][i] == 0
            current[identity] = dict(edge=edge, damage=damage, rest=g['rest_length'][i],
                                     parent=parent, generation=generation,
                                     birth=g['birth_time'][i], stiffness=g['stiffness'][i])
            checks += 1
        assert previous.keys() <= current.keys()
        maximum_stock_error = max(maximum_stock_error, abs(history['stock_error']))
        maximum_energy_fraction = max(maximum_energy_fraction,
                                      abs(history['energy_residual']) / max(1, abs(history['energy']['external_work'])))
        assert abs(history['stock_error']) < 1e-8
        assert abs(history['energy_residual']) < .005 * max(1, abs(history['energy']['external_work']))
        sources.append({'step': frame['step'], 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
        previous = current
    return {'run': str(root), 'run_sha256': hashlib.sha256((root / 'run.json').read_bytes()).hexdigest(),
            'frames': len(sources), 'bond_observations_checked': checks,
            'maximum_stock_error': maximum_stock_error,
            'maximum_energy_residual_over_max_one_work': maximum_energy_fraction,
            'persistent_ids_and_irreversible_damage': True, 'repair_parent_edges_and_generations': True,
            'finite_positions_velocities_and_nonnegative_stock': True, 'source_frames': sources}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('run', type=Path)
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    result = json.dumps(check(a.run), indent=2)
    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(result + '\n')
    print(result)
