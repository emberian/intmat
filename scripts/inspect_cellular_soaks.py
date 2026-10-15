#!/usr/bin/env python3
"""Read-only Linux soak inspection; send over SSH stdin or run on the host.

Reports completed pointers and compact metrics, never loads a Julia checkpoint,
changes a run, or decides that a failed run's last pointer is safe to resume.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def checkpoint_pointer(root):
    for _ in range(3):
        pointer = read_json(root / 'checkpoint.json')
        path = root / pointer['file']
        if path.resolve().parent != root.resolve():
            raise ValueError('checkpoint pointer escapes run directory')
        before = path.stat()
        digest = sha256(path)
        after = path.stat()
        stable = (before.st_ino, before.st_size, before.st_mtime_ns) == (
            after.st_ino, after.st_size, after.st_mtime_ns)
        if stable and pointer == read_json(root / 'checkpoint.json'):
            return dict(pointer=pointer, bytes=after.st_size, sha256=digest,
                        age_seconds=time.time()-after.st_mtime, stable=True)
    return dict(stable=False, note='writer rotated during inspection; retry later')


def metric_chunks(root):
    chunks = sorted(root.glob('metrics-*.jsonl'))
    rows = bad = 0
    recent = []
    last_partial = False
    for chunk in chunks:
        with chunk.open('rb') as stream:
            for raw in stream:
                # A writer may have appended an unfinished last row.
                if not raw.endswith(b'\n'):
                    last_partial = True
                    continue
                try:
                    row = json.loads(raw)
                except (ValueError, UnicodeDecodeError):
                    bad += 1
                    continue
                rows += 1
                m = row.get('metrics', {})
                keep = ('population', 'births', 'deaths', 'bonds', 'blocked_births',
                        'max_generation', 'surviving_lineages', 'biomass',
                        'material_error', 'chemical_energy_error', 'junction_transferred')
                recent.append(dict(step=row.get('step'), time=row.get('time'),
                                   observed_at_unix=row.get('observed_at_unix'),
                                   metrics={k: m[k] for k in keep if k in m}))
                recent = recent[-2:]
    return dict(chunks=len(chunks), complete_rows=rows, invalid_complete_rows=bad,
                partial_tail=last_partial, recent=recent)


def inspect(root):
    root = Path(root)
    now = time.time()
    status = read_json(root / 'status.json')
    result = dict(run=str(root), checked_unix=now, status=status,
                  status_age_seconds=now-status['updated_unix'])
    pid = status['pid']
    proc = Path('/proc') / str(pid)
    result['process_alive'] = proc.is_dir()
    if proc.is_dir():
        result['command'] = (proc / 'cmdline').read_bytes().replace(b'\0', b' ').decode()
        result['cpu_affinity'] = sorted(os.sched_getaffinity(pid))
        result['nice'] = os.getpriority(os.PRIO_PROCESS, pid)
        result['process_state'] = next((x for x in (proc / 'status').read_text().splitlines()
                                       if x.startswith('State:')), None)
        result['process_matches_run'] = (str(root.parent / 'source' / 'scripts' / 'cellular_soak.jl')
                                         in result['command'] and str(root) in result['command'])
        group = next((x[3:] for x in (proc / 'cgroup').read_text().splitlines()
                      if x.startswith('0::')), None)
        if group:
            cg = Path('/sys/fs/cgroup') / group.lstrip('/')
            result['controls'] = {n: (cg/n).read_text().strip() for n in
                                  ('memory.max', 'memory.current', 'memory.swap.max', 'memory.events')
                                  if (cg/n).is_file()}
    result['checkpoint'] = checkpoint_pointer(root)
    result['metrics'] = metric_chunks(root)
    result['artifact_bytes'] = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    source = root.parent / 'source'
    provenance = read_json(source / 'provenance.json')
    mismatches = []
    for name, expected in provenance['hashes'].items():
        path = source / name
        if not path.is_file() or sha256(path) != expected:
            mismatches.append(name)
    result['source_verification'] = dict(files=len(provenance['hashes']), mismatches=mismatches,
                                        manifest_sha256=sha256(source / 'provenance.json'))
    parent = root / 'parent-checkpoint.jls'
    if parent.is_file():
        result['parent_sha256'] = sha256(parent)
    result['failure_files'] = [str(p.relative_to(root)) for p in root.glob('*failure*')]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', nargs='+', help='Existing run directories on this host')
    args = parser.parse_args()
    results = []
    for run in args.runs:
        try:
            results.append(inspect(run))
        except Exception as error:
            results.append(dict(run=run, inspection_error=f'{type(error).__name__}: {error}'))
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
