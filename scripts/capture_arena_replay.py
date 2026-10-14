#!/usr/bin/env python3
"""Capture actual published arena states for a compact, explicitly recorded replay."""
import argparse
import copy
import json
from pathlib import Path
import shutil
import time


def preview(state, side):
    result = copy.deepcopy(state)
    fields = {}
    for name in state['fields']:
        if name in ('vx', 'vy', 'vz', 'activity', 'recovery', 'pressure'):
            continue
        if name not in state['fields']:
            continue
        field = state['fields'][name]
        nx, ny = field['shape']
        stride = max(1, (max(nx, ny) + side - 1) // side)
        xs, ys = list(range(0, nx, stride)), list(range(0, ny, stride))
        fields[name] = {**{key: value for key, value in field.items() if key not in ('shape', 'values')}, 'shape': [len(xs), len(ys)],
                        'values': [round(field['values'][x + nx * y], 4) for y in ys for x in xs]}
    result['fields'] = fields
    for agent in result['agents']:
        for key, value in agent.items():
            if isinstance(value, float):
                agent[key] = round(value, 4)
            elif isinstance(value, list):
                agent[key] = [round(x, 4) if isinstance(x, float) else x for x in value]
    result['preview'] = {'kind': 'recorded actual snapshots', 'field_sampling': 'stride',
                         'maximum_field_side': side, 'display_decimal_places': 4,
                         'note': 'Simulation clock, identities and metrics are retained; display fields and agent floats are sampled or rounded.'}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--frames', type=int, default=25)
    parser.add_argument('--interval', type=float, default=2)
    parser.add_argument('--side', type=int, default=32)
    args = parser.parse_args()
    if args.frames < 1 or args.interval <= 0 or args.side < 1:
        parser.error('frames, interval and side must be positive')
    if args.destination.exists():
        parser.error('choose a fresh destination')
    args.destination.mkdir(parents=True)
    captured = []
    deadline = time.monotonic() + max(30, args.frames * args.interval * 3)
    while len(captured) < args.frames and time.monotonic() < deadline:
        state = json.loads((args.run_directory / 'latest.json').read_text())
        if not captured or state['step'] != captured[-1]['step']:
            captured.append(preview(state, args.side))
            print(f"Captured step {state['step']}, t={state['time']:.3f}", flush=True)
        if len(captured) < args.frames:
            time.sleep(args.interval)
    if not captured:
        raise RuntimeError('no actual snapshots captured')
    for name, data in (('replay.json', captured), ('latest.json', captured[-1])):
        (args.destination / name).write_text(json.dumps(data, separators=(',', ':')))
    shutil.copyfile(Path(__file__).resolve().parents[1] / 'web/arena.html', args.destination / 'arena.html')
    print(f'Saved {len(captured)} actual snapshots to {args.destination}')


if __name__ == '__main__':
    main()
