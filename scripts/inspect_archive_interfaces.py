#!/usr/bin/env python3
"""Measure which designed electrical roots share a saved mineral voxel component.

Six-face connectivity of thresholded cell samples is an observation metric,
not a claim of anatomical connectivity or exact polygon-mesh topology.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np


def inspect(directory, thresholds):
    geometry = json.loads((directory / 'geometry.json').read_text())
    volume = json.loads((directory / 'volume.json').read_text())
    shape = tuple(geometry['shape'])
    assert len(set(shape)) == 1 and volume['order'] == 'F'
    n = shape[0]
    seed = {tuple(x) for x, birth in zip(geometry['positions'], geometry['birth'], strict=True) if birth == 0}
    remaining = set(seed); groups = []
    while remaining:
        queue = [remaining.pop()]; group = set(queue)
        while queue:
            cell = queue.pop()
            for axis in range(3):
                for sign in (-1, 1):
                    other = tuple(x + (sign if d == axis else 0) for d, x in enumerate(cell))
                    if other in remaining:
                        remaining.remove(other); queue.append(other); group.add(other)
        groups.append(group)
    groups.sort(key=min)
    seed_labels = {(x-1) + n*((y-1) + n*(z-1)): label
                   for label, group in enumerate(groups) for x, y, z in group}
    path = directory / 'solid.f32'
    if path.exists():
        solid = np.fromfile(path, dtype='<f4')
    else:
        path = directory / 'solid.f32.gz'
        with gzip.open(path, 'rb') as stream:
            solid = np.frombuffer(stream.read(), dtype='<f4')
    assert solid.size == n**3 and np.isfinite(solid).all()
    results = []
    for threshold in thresholds:
        remaining = set(map(int, np.flatnonzero(solid >= threshold)))
        total = len(remaining); components = []; touches_box = False
        while remaining:
            queue = [remaining.pop()]; count = 0; roots = set()
            while queue:
                index = queue.pop(); count += 1
                if index in seed_labels: roots.add(seed_labels[index])
                coords = (index % n, (index // n) % n, index // (n*n))
                touches_box |= any(x in (0, n-1) for x in coords)
                for delta, coordinate in zip((1, n, n*n), coords):
                    for sign in (-1, 1):
                        if 0 <= coordinate + sign < n:
                            other = index + sign*delta
                            if other in remaining:
                                remaining.remove(other); queue.append(other)
            components.append({'voxels':count, 'seed_root_labels':sorted(roots)})
        results.append({'threshold':threshold, 'foreground_voxels':total,
                        'touches_display_boundary':bool(touches_box),
                        'components':sorted(components, key=lambda x:x['voxels'], reverse=True)})
    return {'event':volume['event'], 'time':volume['time'], 'shape':shape,
            'solid_file_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'seed_root_sizes':[len(group) for group in groups], 'thresholds':results,
            'interpretation':'Six-face connected thresholded mineral samples in the nonperiodic display box. Labels refer to connected components of the original electrical seed. Thermal diffusion itself is periodic; boundary intersection is reported explicitly. This metric is not exact isosurface-mesh topology.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('capture', type=Path)
    parser.add_argument('--thresholds', default='.007,.02,.04')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = inspect(args.capture, [float(x) for x in args.thresholds.split(',')])
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))
