"""Read recorded native fields; occupancy observations and spectra, no simulation."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def native(directory, name, dims):
    return np.fromfile(directory / (name + '.f32'), dtype='<f4').reshape(dims, order='F')


def spectrum(field, extent):
    axes = tuple(i for i, side in enumerate(field.shape) if side > 1)
    shape = tuple(field.shape[i] for i in axes)
    physical = tuple(extent[i] for i in axes)
    field = np.squeeze(field).astype(np.float64)
    fluctuation = field - field.mean()
    power = abs(np.fft.fftn(fluctuation)) ** 2
    frequencies = np.meshgrid(*(2*np.pi*np.fft.fftfreq(n, d=L/n)
                               for n, L in zip(shape, physical)), indexing='ij')
    radius = np.sqrt(sum(k*k for k in frequencies))
    width = 2*np.pi/max(physical)
    shells = np.rint(radius/width).astype(int)
    radial = np.bincount(shells.ravel(), weights=power.ravel())
    radial[0] = 0
    peak = int(np.argmax(radial))
    total = power.sum()
    nyquist = min(np.pi*n/L for n, L in zip(shape, physical))
    return dict(std=float(fluctuation.std()), radial_peak_shell=peak,
                dominant_wavelength=float(2*np.pi/(peak*width)) if peak else None,
                fraction_beyond_point7_axis_nyquist=float(power[radius > .7*nyquist].sum()/total) if total else 0.,
                axis_nyquist=float(nyquist))


def review(root, output, low_fraction=.3, high_fraction=.7):
    root = Path(root); output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    run = json.loads((root/'run.json').read_text())
    dims = tuple(run['dimensions']); extent = tuple(run['extent'])
    capacity = run['parameters']['phase']['capacity']
    low, high = low_fraction*capacity, high_fraction*capacity
    if not 0 < low < high <= capacity:
        raise ValueError('require 0 < low < high <= capacity')
    occupied = np.zeros(dims, bool); ever_opened = np.zeros(dims, bool)
    openings = np.zeros(dims, np.int32); reformations = np.zeros(dims, np.int32)
    selected = set(round((len(run['frame_steps'])-1)*f) for f in (0,.2,.4,.6,.8,1))
    history = []; spectra = []
    for index, step in enumerate(run['frame_steps']):
        directory = root/'frames'/f'{step:07d}'
        b = native(directory, 'biomass', dims); m = native(directory, 'matrix', dims)
        rho = b+m
        opened = occupied & (rho < low)
        filled = (~occupied) & (rho > high)
        reformed = filled & ever_opened
        openings += opened; reformations += reformed
        ever_opened |= opened; occupied[opened] = False; occupied[filled] = True
        row = dict(time=float(step*run['dt']), fraction_below_low=float((rho<low).mean()),
                   fraction_above_high=float((rho>high).mean()),
                   openings=int(openings.sum()), reformations=int(reformations.sum()))
        history.append(row)
        if index in selected:
            spectra.append(dict(time=row['time'], excitation=spectrum(native(directory,'excitation',dims),extent),
                                matrix=spectrum(m,extent), solid=spectrum(rho,extent),
                                biomass_quantiles=np.quantile(b,[.1,.5,.9]).tolist(),
                                matrix_quantiles=np.quantile(m,[.1,.5,.9]).tolist()))
    result = dict(format='living-tissue-native-review-v1', id=run['id'], dimensions=dims, extent=extent,
                  dt=run['dt'], steps=run['steps'], source_sha256=run['source_sha256'],
                  review_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  low=low, high=high,
                  convention='Sampled native rho=B+M. A site must first be observed above high before opening below low; reformation requires an observed opening. Thresholds are concentration fractions of declared solid capacity, not mechanical cavities or organisms. Radial spectra sum fluctuation power in rounded physical wave-number shells; high-frequency fraction uses .7 of the smallest axis Nyquist.',
                  opened_sites=int((openings>0).sum()), reformed_sites=int((reformations>0).sum()),
                  opening_events=int(openings.sum()), reformation_events=int(reformations.sum()),
                  maximum_openings_at_site=int(openings.max()), spectra=spectra, occupancy_history=history,
                  initial=run['history'][0], final=run['history'][-1], elapsed_seconds=run['elapsed_seconds'])
    output.write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root'); parser.add_argument('output')
    parser.add_argument('--low', type=float, default=.3); parser.add_argument('--high', type=float, default=.7)
    args = parser.parse_args()
    result = review(args.root,args.output,args.low,args.high)
    print(json.dumps({k:result[k] for k in ('id','opening_events','reformation_events','opened_sites','reformed_sites')}))
