#!/usr/bin/env python3
"""Compare actual phase simulations on nested index-origin periodic lattices.

No interpolation: refined samples at factor*i align with the old study's
analytic initialization. Central slices retain their actual pixel grids.
"""
import argparse
import array
import hashlib
import json
import math
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('baseline', type=Path)
parser.add_argument('refined', type=Path)
args = parser.parse_args()

def manifest(root):
    return json.loads((root / 'volume' / 'volume.json').read_text())

def phi(root):
    a = array.array('f')
    a.frombytes((root / 'volume' / 'phi.f32').read_bytes())
    return a

def slice_ppm(path, values, n):
    with path.open('wb') as stream:
        stream.write(f'P6\n{n} {n}\n255\n'.encode())
        for value in values:
            u = (1 + math.tanh(2 * value)) / 2
            stream.write(bytes((round(255*u), round(255*math.sqrt(u)), round(255*(1-u)))))

base_shape, high_shape = manifest(args.baseline)['shape'], manifest(args.refined)['shape']
assert len(set(base_shape)) == len(set(high_shape)) == 1
nb, nh = base_shape[0], high_shape[0]
assert nh % nb == 0
factor = nh // nb
base = phi(args.baseline)
assert len(base) == nb**3
squared = largest = 0.
count = positive = 0
finite = True
with (args.refined / 'volume' / 'phi.f32').open('rb') as stream:
    for k in range(nh):
        plane = array.array('f')
        plane.frombytes(stream.read(nh*nh*4))
        assert len(plane) == nh*nh
        finite = finite and all(map(math.isfinite, plane))
        positive += sum(v > 0 for v in plane)
        if k == nh//2:
            slice_ppm(args.refined / 'comparison-refined-slice.ppm', plane, nh)
        if k % factor:
            continue
        kb = k//factor
        for jb in range(nb):
            for ib in range(nb):
                difference = plane[ib*factor + jb*factor*nh] - base[ib + jb*nb + kb*nb*nb]
                squared += difference*difference
                largest = max(largest, abs(difference))
                count += 1
    assert stream.read(1) == b''
slice_ppm(args.refined / 'comparison-baseline-slice.ppm', base[(nb//2)*nb*nb:(nb//2+1)*nb*nb], nb)
metadata = [json.loads((root/'run.json').read_text()) for root in (args.baseline,args.refined)]
assert metadata[1]['extent'] == [64,64,64]
assert metadata[0]['traces'][-1]['extent'] == [64,64,64]
assert metadata[0]['traces'][-1]['history'][-1]['time'] == metadata[1]['time']
assert metadata[0]['traces'][-1]['dt'] == metadata[1]['dt']
result = dict(baseline=str(args.baseline), refined=str(args.refined),
              shapes=[base_shape,high_shape], spatial_factor_per_axis=factor,
              cell_count_factor=factor**3, physical_extent=[64,64,64],
              nominal_cell_spacing=[64/nb,64/nh], time=metadata[1]['time'],
              comparison='Matching index-origin lattice samples, not interpolated or block-averaged. Both simulations independently initialized and evolved.',
              sampled_rms_difference=math.sqrt(squared/count), sampled_max_difference=largest,
              refined_phi_finite=finite, refined_positive_fraction=positive/nh**3,
              run_sha256=[hashlib.sha256((root/'run.json').read_bytes()).hexdigest() for root in (args.baseline,args.refined)])
result['discrete_energy'] = [metadata[0]['traces'][-1]['history'][-1]['energy'], metadata[1]['observations'][-1]['energy']]
result['relative_discrete_energy_difference'] = result['discrete_energy'][1]/result['discrete_energy'][0]-1
result['all_exported_raw_lengths_match_native_shape'] = all(
    (args.refined/'volume'/f'{name}.f32').stat().st_size == 4*nh**3
    for name in manifest(args.refined)['fields'])
assert result['all_exported_raw_lengths_match_native_shape']
(args.refined/'resolution-comparison.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
