#!/usr/bin/env python3
"""Package bounded actual per-cell program captures; never reconstruct outcomes."""
import argparse, hashlib, json, math, pathlib, re

SCHEMA = 'interior-cell-program-observatory-1'

def validate(document):
    if document.get('schema') != SCHEMA:
        raise ValueError('Expected actual cell-program observatory schema')
    frames = document.get('frames')
    if not isinstance(frames, list) or not frames:
        raise ValueError('At least one completed capture is required')
    last_time, last_step = -1, -1
    for frame in frames:
        time, step, cells = frame.get('time'), frame.get('step'), frame.get('cells')
        if type(time) not in (int, float) or not math.isfinite(time) or time < last_time or type(step) is not int or step < last_step:
            raise ValueError('Actual capture clocks must be finite and nondecreasing')
        if not isinstance(cells, list) or len(cells) > 256:
            raise ValueError('Bounded actual resident records required')
        ids = [c.get('id') for c in cells]
        if any(type(i) is not int or i <= 0 for i in ids) or len(set(ids)) != len(ids):
            raise ValueError('Unique positive resident IDs required per capture')
        for cell in cells:
            for key in ('raw', 'food', 'catalyst', 'fast_raw', 'slow_raw', 'raw_next', 'prediction_error', 'funded', 'requested'):
                values = cell.get(key)
                if values is not None and (not isinstance(values, list) or len(values) != 2 or any(type(x) not in (int, float) or not math.isfinite(x) for x in values)):
                    raise ValueError(f'{key} must contain two actual finite values or be absent/null')
            if 'position' in cell and (not isinstance(cell['position'], list) or len(cell['position']) != 3 or any(type(x) not in (int, float) or not math.isfinite(x) for x in cell['position'])):
                raise ValueError('Position requires three actual finite coordinates')
        last_time, last_step = time, step
    return document

def build(recordings, destination, max_bytes=4_000_000, run_report=None):
    runs = []
    for source in recordings:
        raw = source.read_bytes()
        document = validate(json.loads(raw))
        runs.append({'recording_sha256': hashlib.sha256(raw).hexdigest(), **document})
    by_arm = {r.get('metadata', {}).get('arm'): r for r in runs}
    note = None
    if 'predictive' in by_arm and 'persistence' in by_arm:
        keys = ('id', 'program', 'phase', 'investment', 'allocation', 'reserve', 'funded', 'requested', 'position', 'mass', 'energy', 'stores')
        def observed_choices(run):
            return [(f['step'], f['time'], [{k: c.get(k) for k in keys} for c in f['cells']]) for f in run['frames']]
        if observed_choices(by_arm['predictive']) == observed_choices(by_arm['persistence']):
            note = 'Authored controller baseline. Predictive and current-reading controllers chose the same actions.'
            if run_report is not None:
                report = json.loads(run_report.read_text())
                if len(runs) == 3 and all(abs(report[r['metadata']['arm']]['final']['time'] - 8) < 1e-8 for r in runs):
                    note = 'Authored controller baseline · three eight-unit runs. Predictive and current-reading controllers chose the same actions.'
    payload = json.dumps({'schema': 'cell-program-explorer-1', 'comparison_note': note, 'runs': runs}, separators=(',', ':'), allow_nan=False)
    if re.search(r'/Users/|/home/|/othersys/|/tank/|/private/|/tmp/|file://|localhost:|127\.0\.0\.1:', payload):
        raise ValueError('Public recording includes private paths or loopback locations')
    template = pathlib.Path(__file__).with_name('cell_program_observatory.html').read_text()
    html = template.replace('__DATA__', payload.replace('</', '<\\/'))
    if len(html.encode()) > max_bytes:
        raise ValueError('Viewer exceeds declared byte budget; select fewer IDs/captures explicitly')
    destination.mkdir(parents=True, exist_ok=True)
    (destination / 'index.html').write_text(html)
    return {'bytes': len(html.encode()), 'runs': len(runs), 'captures': sum(len(r['frames']) for r in runs)}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--recording', action='append', type=pathlib.Path, required=True)
    p.add_argument('--destination', type=pathlib.Path, required=True)
    p.add_argument('--run-report', type=pathlib.Path, help='Optional measured endpoint report for run-duration caption')
    p.add_argument('--max-bytes', type=int, default=4_000_000)
    a = p.parse_args()
    print(json.dumps(build(a.recording, a.destination, a.max_bytes, a.run_report)))

if __name__ == '__main__': main()
