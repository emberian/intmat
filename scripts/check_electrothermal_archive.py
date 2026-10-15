#!/usr/bin/env python3
"""Check saved electrical identities, face topology and final material balance."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def validate(root):
    run_path = root / 'run.json'
    run = json.loads(run_path.read_text())
    shape = tuple(run['shape'])
    previous = {}
    observations = 0
    records = []
    initial_count = None
    for frame, history in zip(run['frames'], run['history'], strict=True):
        directory = root / frame['directory']
        volume = json.loads((directory / 'volume.json').read_text())
        assert tuple(volume['shape']) == shape and volume['order'] == 'F' and volume['dtype'] == 'float32'
        assert frame['step'] == history['event'] == volume['event']
        assert abs(frame['time'] - volume['time']) <= 1e-5 * max(1, frame['time'])
        if records:
            assert frame['step'] > records[-1]['step'], 'Non-increasing captured events'
        gpath = directory / 'geometry.json'
        g = json.loads(gpath.read_text())
        positions = [tuple(p) for p in g['positions']]
        assert len(set(positions)) == len(positions), 'Duplicate occupied cell'
        birth = dict(zip(positions, g['birth'], strict=True))
        assert all(len(p) == 3 and all(1 <= p[d] <= shape[d] for d in range(3)) for p in positions)
        assert tuple(g['shape']) == shape
        assert all(p in birth and birth[p] == b for p, b in previous.items()), 'Changed historical conductor'
        nonzero_births = sorted(b for b in birth.values() if b)
        assert nonzero_births == list(range(1, frame['step'] + 1)), 'Each event must add exactly one cell'
        initial_count = initial_count if initial_count is not None else len(positions)
        assert len(positions) == initial_count + frame['step'] == history['occupied']
        index = {p: i + 1 for i, p in enumerate(positions)}
        expected = set()
        for p, i in index.items():
            for axis in range(3):
                q = tuple(x + int(d == axis) for d, x in enumerate(p))
                if q in index:
                    expected.add(tuple(sorted((i, index[q]))))
        actual = {tuple(sorted(e)) for e in g['edges']}
        assert len(actual) == len(g['edges']) and actual == expected, 'Saved edges differ from exact face adjacency'
        neighbors = {p: [] for p in positions}
        for i, j in actual:
            neighbors[positions[i-1]].append(positions[j-1]); neighbors[positions[j-1]].append(positions[i-1])
        assert all(b == 0 or any(birth[q] < b for q in neighbors[p]) for p, b in birth.items()), 'Growth lacks an older face neighbor'
        remaining = set(positions); sizes = []
        while remaining:
            queue = [remaining.pop()]; count = 0
            while queue:
                p = queue.pop(); count += 1
                for q in neighbors[p]:
                    if q in remaining:
                        remaining.remove(q); queue.append(q)
            sizes.append(count)
        assert history['finite'] and history['nutrient_min'] >= 0
        assert abs(history['material_relative_error']) < 1e-6
        assert history['residual'] <= run['params']['tolerance']
        records.append({'step':frame['step'], 'occupied':len(positions), 'components':len(sizes),
                        'component_sizes':sorted(sizes, reverse=True), 'geometry_sha256':digest(gpath)})
        previous = birth; observations += len(positions)
    assert run['frames'][-1]['step'] == run['steps']
    final = root / run['frames'][-1]['directory']
    fields = {}
    field_hashes = {}
    for name in ('temperature', 'dose', 'solid'):
        path = final / (name + '.f32')
        if path.exists():
            field = np.fromfile(path, dtype='<f4')
        else:
            path = final / (name + '.f32.gz')
            with gzip.open(path, 'rb') as stream:
                field = np.frombuffer(stream.read(), dtype='<f4')
        assert field.size == np.prod(shape) and np.isfinite(field).all() and field.min() >= 0
        fields[name] = field
        field_hashes[name] = digest(path)
    food_path = root / 'final-nutrient.f32'
    food = np.fromfile(food_path, dtype='<f4')
    assert food.size == np.prod(shape) and np.isfinite(food).all() and food.min() >= 0
    local_error = float(np.max(np.abs(food.astype(np.float64) + fields['solid'] - 1)))
    mass_error = float((food.sum(dtype=np.float64) + fields['solid'].sum(dtype=np.float64)) / food.size - 1)
    assert local_error < 5e-5 and abs(mass_error) < 1e-6
    return {'status':'passed', 'run_sha256':digest(run_path), 'shape':shape, 'steps':run['steps'],
            'termination':run.get('termination','completed'), 'captures':len(records),
            'conductor_observations':observations, 'initial_components':records[0]['components'],
            'final_components':records[-1]['components'], 'final_pointwise_material_error':local_error,
            'final_relative_material_error':mass_error, 'final_field_sha256':field_hashes,
            'final_nutrient_sha256':digest(food_path), 'frames':records}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('run', type=Path); p.add_argument('--output', type=Path)
    a = p.parse_args(); report = validate(a.run)
    output = a.output or a.run / 'validation.json'
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'frames'}))
