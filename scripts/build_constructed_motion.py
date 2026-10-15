#!/usr/bin/env python3
"""Pack actual 2D constructive observations into a small, standalone motion exhibit."""
import argparse
import array
import base64
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
FIELDS = ['enzyme_a', 'enzyme_b', 'retained_a', 'retained_b']
SCALARS = ['energy', 'energy_capacity', 'biomass']

def packed(kind, values):
    data = array.array(kind, values)
    if sys.byteorder != 'little': data.byteswap()
    return base64.b64encode(data.tobytes()).decode()

def quant(value, maximum, levels):
    return max(0, min(levels, round(float(value) / (maximum or 1) * levels)))

def build(source, destination):
    source, destination = Path(source), Path(destination)
    runs = []
    scales = {key: 0. for key in FIELDS + SCALARS + ['radius']}
    # Two streaming passes avoid retaining hundreds of large source dictionaries.
    for seed in [211, 307]:
        directory = source / f'strong-{seed}-continuation-3100'
        files = sorted(directory.glob('frame-*.json'))
        if len(files) != 201: raise ValueError(f'{seed}: expected 201 actual frames')
        identities = {}
        for path in files:
            frame = json.loads(path.read_text())
            if frame['grid']['shape'][2] != 1: raise ValueError('This exhibit is explicitly 2D')
            for key in FIELDS: scales[key] = max(scales[key], max(frame['fields'][key]['values']))
            for cell in frame['agents']:
                identities[cell['id']] = [cell['id'], cell['parent'], cell['lineage'], cell['generation']]
                for key in SCALARS + ['radius']: scales[key] = max(scales[key], cell[key])
        runs.append((seed, directory, files, sorted(identities.values())))
    result = {'schema': 'constructed-motion-1', 'endian': 'little', 'scales': scales,
              'fields': FIELDS, 'field_encoding': 'temporal-delta-u8-mod256', 'scalar_order': SCALARS, 'runs': []}
    provenance = {'title': 'Kept, released, remade.', 'family': 'ecology',
                  'subtitle': 'Two recordings · actual 2D states',
                  'description': 'Catalyst and retained food surround physical cells with persistent bonds. Watch two completed constructive ecologies, select any recorded cell, and compare its neural state with paid material and work flows.',
                  'viewer': 'index.html', 'viewer_label': 'Open recorded motion',
                  'viewer_context': 'Recorded motion + cell inspection',
                  'image_caption': 'Actual planar state · seed 211 · t 3064',
                  'native_shape': [32, 32, 1], 'kind': 'compact display of completed 2D recordings and exact-replay paid-flow sidecars; no altered trajectories or temporal cell interpolation',
                  'encoding': {'fields': 'unsigned 8-bit linear, per-field maxima shared across both recordings; lossless temporal byte differences',
                               'positions_radius_energy_capacity_biomass': 'unsigned 16-bit linear',
                               'neural': 'signed 8-bit linear [-1,1], max error 1/254',
                               'metabolites_and_paid_flows': 'IEEE Float32, rounded from recorded Float64'},
                  'runs': []}
    for seed, directory, files, identities in runs:
        index = {row[0]: i for i, row in enumerate(identities)}
        flowdir = source / 'analysis' / f'paid-flow-{seed}'
        report = json.loads((flowdir / 'report.json').read_text())
        if not report['complete']: raise ValueError(f'{seed}: paid-flow replay incomplete')
        flowfile = flowdir / 'interval-cell-flows.jsonl'
        flows = {row['step']: row for row in map(json.loads, flowfile.read_text().splitlines())}
        run = {'seed': seed, 'identities': identities, 'frames': [],
               'focus': {'ids': [728, 767] if seed == 211 else [917, 983],
                         'window': [3061, 3069] if seed == 211 else [3081.5, 3089.5]}}
        previous_fields = {}
        for path in files:
            raw = path.read_bytes(); frame = json.loads(raw); flow = flows[frame['step']]
            if not flow['physical_snapshot_exact_match']: raise ValueError(f'{path.name}: replay mismatch')
            if hashlib.sha256(raw).hexdigest() != flow['recorded_frame_sha256']: raise ValueError('Frame identity mismatch')
            if len(frame['agents']) != frame['metrics']['population']: raise ValueError('Population was truncated')
            run['grid'] = frame['grid']
            local = {row['id']: row for row in flow['cells']}
            agents = frame['agents']; n = len(agents)
            geometry = []; scalar = []; stores = []; neural = []; allocation = []; paid = []
            for cell in agents:
                if not (len(cell['hidden']) == len(cell['adaptation']) == 8 and len(cell['actions']) == 11): raise ValueError('Unsupported controller shape')
                pos = cell['cell_position']; shape = frame['grid']['shape']
                geometry.extend([quant(pos[0]/shape[0], 1, 65535), quant(pos[1]/shape[1], 1, 65535), quant(cell['radius'], scales['radius'], 65535)])
                scalar.extend(quant(cell[key], scales[key], 65535) for key in SCALARS)
                stores.extend(cell['metabolites'])
                neural.extend(max(-127, min(127, round(v*127))) for key in ['hidden', 'adaptation', 'actions'] for v in cell[key])
                allocation.extend(quant(v, 1, 255) for v in cell['metabolism']['allocation'] + [cell['catalysis']['investment']] + cell['catalysis']['allocation'])
                f = local.get(cell['id'])
                paid.extend((f['actual_secretion'] + [f['actual_construction_work'], f['actual_motor_work']]) if f else [float('nan')]*4)
            fields = {}
            for key in FIELDS:
                current = [quant(v, scales[key], 255) for v in frame['fields'][key]['values']]
                previous = previous_fields.get(key, [0] * len(current))
                fields[key] = packed('B', ((a-b) % 256 for a, b in zip(current, previous)))
                previous_fields[key] = current
            interval = frame['observation_interval']
            run['frames'].append({'time': frame['time'], 'ids': packed('H', (index[c['id']] for c in agents)),
                                  'geometry': packed('H', geometry), 'scalar': packed('H', scalar),
                                  'stores': packed('f', stores), 'neural': packed('b', neural),
                                  'allocation': packed('B', allocation), 'paid': packed('f', paid),
                                  'bonds': packed('H', (index[i] for edge in frame['metrics']['cell_bonds'] for i in edge)),
                                  'fields': fields,
                                  'interval': {'duration': interval['duration'], 'total_secretion': flow['actual_total_secretion'],
                                               'unassigned': flow['secretion_not_attributed_to_observed_survivors'],
                                               'food_solubilized': interval['actual_flows']['food_solubilized']}})
        result['runs'].append(run)
        sequence = directory / 'sequence.json'
        provenance['runs'].append({'seed': seed, 'frames': len(files), 'start': run['frames'][0]['time'], 'end': run['frames'][-1]['time'],
                                   'sequence_sha256': hashlib.sha256(sequence.read_bytes()).hexdigest(),
                                   'paid_flow_sha256': hashlib.sha256(flowfile.read_bytes()).hexdigest(),
                                   'paid_flow_report_sha256': hashlib.sha256((flowdir/'report.json').read_bytes()).hexdigest(),
                                   'all_physical_snapshots_exact_match': True, 'focus': run['focus']})
    destination.mkdir(parents=True, exist_ok=True)
    data = gzip.compress(json.dumps(result, separators=(',', ':'), allow_nan=False).encode(), compresslevel=9, mtime=0)
    (destination / 'recording.json.gz').write_bytes(data)
    provenance['recording_sha256'] = hashlib.sha256(data).hexdigest()
    template = ROOT / 'web' / 'constructed_motion.html'
    shutil.copyfile(template, destination / 'index.html')
    provenance['viewer_sha256'] = hashlib.sha256(template.read_bytes()).hexdigest()
    (destination / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    size = sum(path.stat().st_size for path in destination.iterdir() if path.is_file())
    print(json.dumps({'destination': str(destination), 'recording_bytes': len(data), 'total_bytes': size}))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT / 'outputs/constructed-feedback-observation-20261015')
    parser.add_argument('--destination', type=Path, default=ROOT / 'outputs/constructed-motion-20261015')
    args = parser.parse_args(); build(args.source, args.destination)
