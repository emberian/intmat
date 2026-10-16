#!/usr/bin/env python3
"""Package recorded continuous-ramp probes; no simulation or inferred recordings."""
import argparse, hashlib, json, math, pathlib, re, statistics

FILES = {'authority': 'continuous-ramp-authority.json', 'corrected9301': 'continuous-ramp-construction.json',
         'corrected9403': 'continuous-ramp-construction-noise9403.json', 'saturated9301': 'continuous-ramp-saturated-critic.json'}

def pack(source):
    reports = {k: json.loads((source / name).read_text()) for k, name in FILES.items()}
    out = {'schema': 'continuous-ramp-explorer-1', 'authority': [], 'runs': {},
           'source_hashes': {name: hashlib.sha256((source / name).read_bytes()).hexdigest() for name in FILES.values()}}
    keys = ('reward', 'prediction', 'drive', 'funded', 'balance', 'trend', 'gain')
    for row in reports['authority']['results']:
        assert len(row['reward']) == 4 and row['balance'] == [0.] * 4
        out['authority'].append({k: row[k] for k in ('direction', 'logit', 'supply_duration', 'reward', 'funded', 'supplied', 'trend')})
    assert len(out['authority']) == 12
    for name, report in reports.items():
        if name == 'authority': continue
        run = {k: report[k] for k in ('seed', 'noise_seed', 'episodes', 'critic_scale', 'validation')}
        assert report['episodes'] == 64 and report['validation']['funded_spread'] < 1e-7
        run['arms'] = {}
        for arm, value in report['arms'].items():
            probes = []
            assert len(value['probes']) == 8
            for row in value['probes']:
                assert len(row['reward']) == len(row['drive']) == 4
                if 'prediction' in row: assert len(row['prediction']) == 4
                assert all(abs(v) < 1e-7 for v in row['balance'])
                if row['context'] <= 2: assert row['balance'] == [0.] * 4
                assert all(abs(sum(row['features'][i + 4*k]**2 for k in range(3)) - 2) < 1e-6 for i in range(4))
                probes.append({k: row[k] for k in ('context', 'swap', *keys) if k in row})
                baseline = next(p for p in report['arms']['no_credit']['probes'] if p['context'] == row['context'] and not p['swap'])
                assert abs(row['gain'] - statistics.mean(row['reward']) + statistics.mean(baseline['reward'])) < 1e-12
            run['arms'][arm] = {'probes': probes, 'learned': value['learned'], 'work': value['total_learning_work']}
        out['runs'][name] = run
    payload = json.dumps(out, separators=(',', ':'), allow_nan=False)
    assert not re.search(r'/Users/|/home/|file://|localhost:|127\.0\.0\.1:', payload)
    return payload

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=pathlib.Path, default=pathlib.Path('docs/figures'))
    parser.add_argument('--destination', type=pathlib.Path, default=pathlib.Path('site/laboratory/continuous-ramp'))
    args = parser.parse_args()
    payload = pack(args.source)
    html = pathlib.Path(__file__).with_name('ramp_construction_viewer.html').read_text().replace('__DATA__', payload.replace('</', '<\\/'))
    assert len(html.encode()) < 90_000, 'Keep this explanation a small static instrument'
    args.destination.mkdir(parents=True, exist_ok=True)
    (args.destination / 'index.html').write_text(html)
    print(f'{len(html.encode())} HTML bytes, including {len(payload)} measured-data bytes')

if __name__ == '__main__': main()
