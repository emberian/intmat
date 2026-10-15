#!/usr/bin/env python3
"""Summarize retained arena replay or metabolic study JSON; never infer lost history."""
import argparse
import json
from pathlib import Path


def summarize(path):
    data = json.loads(Path(path).read_text())
    frames = data if isinstance(data, list) else data.get('samples', [data])
    def metrics(frame):
        return frame.get('metrics', frame)
    alive = [f for f in frames if metrics(f).get('population', 0) > 0]
    extinct = next((f for f in frames if metrics(f).get('population') == 0), None)
    result = dict(path=str(path), retained_samples=len(frames),
                  time_range=[frames[0]['time'], frames[-1]['time']],
                  max_sampled_population=max(metrics(f).get('population', 0) for f in frames),
                  first_sampled_extinction=extinct['time'] if extinct else None)
    # These are bounds on retained observations, not whole-run extrema.
    if alive:
        result['last_alive_time'] = alive[-1]['time']
        result['last_alive_metrics'] = {k: v for k, v in metrics(alive[-1]).items()
                                       if k not in ('morphology',)}
    result['last_metrics'] = metrics(frames[-1])
    if data and isinstance(data, dict):
        result['condition'] = data.get('condition')
        result['seed'] = data.get('seed')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='+')
    args = parser.parse_args()
    print(json.dumps([summarize(p) for p in args.paths], indent=2))


if __name__ == '__main__':
    main()
