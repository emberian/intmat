"""Resolve the actual inputs read by collection renderers for cache validation."""
import hashlib
import json
from pathlib import Path


SAVED_FRAME_RENDERERS = {
    'interfacial-choir': 'render_interfacial_choir.py',
    'electrothermal-archive': 'render_electrothermal_archive.py',
    'repair-weave': 'render_repair_weave.py',
    'constructive-volume': 'render_constructive_volume.py',
}


def saved_frame_inputs(directory, job):
    """Follow each renderer's saved-frame selection and raw/gzip precedence."""
    directory = Path(directory)
    manifest = directory / 'run.json'
    run = json.loads(manifest.read_text())
    args = job.get('args', [])
    selection = 'all'
    for index, argument in enumerate(args):
        if argument == '--steps':
            selection = args[index + 1]
        elif argument.startswith('--steps='):
            selection = argument.partition('=')[2]
    steps = None if selection == 'all' else set(map(int, selection.split(',')))
    frames = [frame for frame in run['frames'] if steps is None or frame['step'] in steps]
    if not frames:
        raise ValueError('No saved frames match the render selection')
    inputs = [manifest]
    for frame in frames:
        source = directory / frame['directory']
        if job.get('renderer') == 'repair-weave':
            inputs.append(source / 'geometry.json')
        elif job.get('renderer') == 'constructive-volume':
            inputs.extend(source / name for name in ('volume.json', 'geometry.json', 'retained_a.f32', 'retained_b.f32', 'permeability.f32'))
        elif job.get('renderer') == 'electrothermal-archive':
            inputs.extend(source / name for name in ('volume.json', 'geometry.json'))
            for name in ('solid.f32', 'temperature.f32'):
                raw = source / name
                inputs.append(raw if raw.exists() else source / (name + '.gz'))
        else:
            inputs.extend(source / name for name in ('volume.json', 'solid.f32', 'activity.f32', 'pressure.f32'))
    return inputs


def source_revision(directory, job):
    """Cheap watch key; ctime also detects rewrites with restored mtime."""
    directory = Path(directory)
    inputs = saved_frame_inputs(directory, job) if job.get('renderer') in SAVED_FRAME_RENDERERS else [
        path for path in sorted((directory / job.get('volume', 'volume')).glob('*'))
        if path.suffix in ('.f32', '.json') and path.name != 'optical-mapping.json']
    return tuple((str(path), path.stat().st_size, path.stat().st_mtime_ns,
                  path.stat().st_ctime_ns) for path in inputs)


def source_signature(directory, job):
    directory = Path(directory)
    if job.get('renderer') in SAVED_FRAME_RENDERERS:
        records = []
        for path in saved_frame_inputs(directory, job):
            digest = hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
                    digest.update(block)
            records.append((str(path.relative_to(directory)), path.stat().st_size, digest.hexdigest()))
        return records
    volume = directory / job.get('volume', 'volume')
    return [(path.name, path.stat().st_size, path.stat().st_mtime_ns)
            for path in sorted(volume.glob('*'))
            if path.suffix in ('.f32', '.json') and path.name != 'optical-mapping.json']


# Retain the previous helper name for external queue consumers.
interfacial_inputs = saved_frame_inputs
