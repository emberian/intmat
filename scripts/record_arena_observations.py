#!/usr/bin/env python3
"""Read-only HTTP recorder for published arena states; never reconstructs missed history."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import json
import math
import os
from pathlib import Path
import signal
import time
import urllib.request
import uuid


def encoded(value):
    return (json.dumps(value, separators=(',', ':'), allow_nan=False) + '\n').encode()


def summarize(state):
    agents = state['agents']
    signatures = Counter(a.get('genome_signature') for a in agents)
    n = len(agents)
    traits = {}
    for key in ('energy', 'storage', 'plastic_norm', 'generation', 'age'):
        values = [a[key] for a in agents if isinstance(a.get(key), (int, float)) and math.isfinite(a[key])]
        if values:
            traits[key] = {'n': len(values), 'mean': sum(values) / len(values), 'min': min(values), 'max': max(values)}
    return {'step': state['step'], 'model_time': state['time'], 'generated_at_unix': state.get('generated_at_unix'),
            'metrics': state['metrics'], 'captured_agents': n, 'actual_population': state['metrics']['population'],
            'agents_complete': n == state['metrics']['population'], 'display_cap': 3000,
            'display_cap_source': 'src/arena.jl arena_snapshot default; HTTP schema does not advertise cap',
            'captured_traits': traits, 'captured_signature_counts': dict(signatures),
            'captured_signature_diversity': len(signatures),
            'captured_signature_entropy': -sum((v/n)*math.log(v/n) for v in signatures.values()) if n else 0,
            'paused': state.get('paused'), 'stopped': state.get('stopped')}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', action='append', required=True, help='NAME=HTTP_URL (repeatable)')
    p.add_argument('--output', type=Path, required=True, help='Fresh directory; existing directories are refused')
    p.add_argument('--interval', type=float, default=20)
    p.add_argument('--archive-interval', type=float, default=300)
    p.add_argument('--budget-bytes', type=int, default=6_000_000_000)
    p.add_argument('--timeout', type=float, default=10)
    p.add_argument('--polls', type=int, default=0, help='0 means continuous')
    args = p.parse_args()
    if min(args.interval, args.archive_interval, args.budget_bytes, args.timeout) <= 0 or args.polls < 0:
        p.error('Intervals, budget and timeout must be positive; polls must be nonnegative')
    sources = dict(s.split('=', 1) for s in args.source)
    if len(sources) != len(args.source) or any(not n.replace('-', '').replace('_', '').isalnum() or not u.startswith(('http://', 'https://')) for n,u in sources.items()):
        p.error('Use unique simple source names and HTTP(S) URLs')
    args.output.mkdir(parents=True, exist_ok=False)
    stop = False
    def halt(signum, frame):
        nonlocal stop
        stop = True
    signal.signal(signal.SIGTERM, halt)
    signal.signal(signal.SIGINT, halt)
    status = {'pid': os.getpid(), 'started_at': time.time(), 'state': 'running', 'sources': {}, 'bytes_written': 0}
    def atomic_status():
        status['updated_at'] = time.time()
        tmp = args.output / 'status.tmp'
        tmp.write_bytes(encoded(status))
        tmp.replace(args.output / 'status.json')
    def write(path, data, append=False):
        # Reserve space for status and errors, and never remove old observations.
        if status['bytes_written'] + len(data) > args.budget_bytes - 65536:
            raise OverflowError('Storage budget reached; recording stopped without deleting data')
        with path.open('ab' if append else 'xb') as f:
            f.write(data)
        status['bytes_written'] += len(data)
    manifest = {'pid': os.getpid(), 'parameters': {**vars(args), 'output': str(args.output.resolve())}, 'sources': sources,
                'started_at_unix': status['started_at'], 'scope': 'Published HTTP observations only; no reconstruction of earlier history'}
    write(args.output / 'manifest.json', encoded(manifest))
    previous, archived = {}, {}
    polls = 0
    try:
        while not stop and (not args.polls or polls < args.polls):
            for name, url in sources.items():
                if stop:
                    break
                now = time.time()
                try:
                    with urllib.request.urlopen(url, timeout=args.timeout) as response:
                        raw = response.read(64 * 1024 * 1024 + 1)
                    if len(raw) > 64 * 1024 * 1024:
                        raise ValueError('Snapshot exceeds 64 MiB read limit')
                    state = json.loads(raw)
                    row = summarize(state)
                    identity = (state['step'], state.get('generated_at_unix'))
                    current = status['sources'].setdefault(name, {})
                    current.update(last_success_at=now, step=state['step'], generated_at_unix=state.get('generated_at_unix'),
                                   snapshot_age_seconds=now-state['generated_at_unix'] if state.get('generated_at_unix') else None)
                    if previous.get(name) != state['step']:
                        row.update(observed_at_unix=now, observed_at_utc=datetime.fromtimestamp(now, timezone.utc).isoformat(), source=name)
                        write(args.output / (name + '.jsonl'), encoded(row), append=True)
                        previous[name] = state['step']
                        current['last_distinct_step_at'] = now
                    if identity != archived.get(name, (None,))[0] and now - archived.get(name, (None, 0))[1] >= args.archive_interval:
                        filename = f'{name}-{state["step"]}-{time.time_ns()}-{uuid.uuid4().hex[:8]}.json.gz'
                        write(args.output / filename, gzip.compress(raw, mtime=0))
                        archived[name] = (identity, now)
                        current['latest_archive'] = filename
                    current['last_error'] = None
                except OverflowError:
                    raise
                except Exception as exc:
                    error = {'observed_at_unix': now, 'source': name, 'error': f'{type(exc).__name__}: {exc}'}
                    write(args.output / 'errors.jsonl', encoded(error), append=True)
                    status['sources'].setdefault(name, {})['last_error'] = error
                atomic_status()
            polls += 1
            if not stop and (not args.polls or polls < args.polls):
                until = time.monotonic() + args.interval
                while not stop and time.monotonic() < until:
                    time.sleep(min(.25, max(0, until-time.monotonic())))
        status['state'] = 'stopped' if stop else 'completed'
    except BaseException as exc:
        status.update(state='budget_exhausted' if isinstance(exc, OverflowError) else 'failed', error=f'{type(exc).__name__}: {exc}')
        if not isinstance(exc, OverflowError):
            raise
    finally:
        atomic_status()


if __name__ == '__main__':
    main()
