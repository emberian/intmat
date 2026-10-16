"""Prepare actual-state frames for the shared animation curator; never encode here.

Writes sequence.json and a pending review candidate beside the rendered PNGs.
After visual review, use scripts/curate_animation.py for encoding/publication.
Existing sequence/review files are deliberately not overwritten.
"""
import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('capture', type=Path)
    p.add_argument('frames', type=Path)
    p.add_argument('--id', default='voltage-pieces-motion-48')
    p.add_argument('--title', default='The field gives way')
    p.add_argument('--fps', type=float, default=8)
    a = p.parse_args()
    capture = json.loads((a.capture / 'capture.json').read_text())
    record = json.loads((a.frames / 'frames.json').read_text())
    verify = json.loads((a.capture / 'verification.json').read_text())
    sequence_path = a.frames / 'sequence.json'
    review_path = a.frames / 'art-reviewed.candidate.json'
    if sequence_path.exists() or review_path.exists():
        raise FileExistsError('Preserve existing sequence/review; use a fresh frame directory')
    indices = sorted(map(int, record['frames']))
    if len(indices) < 2 or not 0 < a.fps <= 60:
        raise ValueError('At least two actual frames and fps in (0,60] required')
    frames = []
    for i in indices:
        r, s = record['frames'][str(i)], capture['frames'][i]
        assert (r['step'], r['time'], r['source_sha256']) == (s['step'], s['time'], s['sha256'])
        assert digest(a.frames / r['image']) == r['image_sha256']
        frames.append({'file': r['image'], 'time': r['time'], 'sha256': r['image_sha256']})
    assert all(b['time'] > a['time'] for a, b in zip(frames, frames[1:]))
    base = json.loads((a.capture / 'base.json').read_text())
    source = {
        'study': 'studies/voltage_between_pieces.jl',
        'render_study': 'scripts/render_voltage_pieces.py',
        'actual_samples': len(frames), 'raw_captured_samples': len(capture['frames']),
        'raw_indices': indices, 'native_shape': base['dims'],
        'native_nodes': capture['native_nodes'], 'native_bonds': capture['native_bonds'],
        'temporal_interpolation': False, 'restart_checkpoints': False,
        'camera': 'fixed throughout: (2.8,-4,2.6), ortho2.85', 'samples': record['samples'],
        'potential_levels': [.12, .22, .32, .42, .52, .62, .72, .82, .90],
        'global_dose_optical_max': capture['global_dose_optical_max'],
        'physical_time_interval': [frames[0]['time'], frames[-1]['time']],
        'capture_verification_sha256': digest(a.capture / 'verification.json'),
        'render_record_sha256': digest(a.frames / 'frames.json'),
        'renderer_sha256': record['renderer_sha256'],
        'simulation_source_sha256': capture['source_sha256'],
        'base_sha256': capture['base_sha256'],
        'numerical_closure': {k: verify[k] for k in (
            'maximum_independent_residual', 'maximum_relative_power_defect',
            'maximum_mechanical_relative_work_residual', 'irreversible_damage_and_nondecreasing_dose')},
        'interpretation': 'Authored grains, quenched geometry, weak seams and grip forcing; evolved spring dynamics, irreversible failure, voltage redistribution and Joule dose. Gold foils observe potential iso-surfaces, not extra physical sheets. Fixed optical transfer/camera. No temporal interpolation or closed thermodynamic claim.',
    }
    sequence = {'id': a.id, 'title': a.title,
                'description': 'A fractured support opens and redirects its electrical interior. Gold potential contours separate as color records accumulated Joule dose. Actual Julia observations under one fixed camera and optical transfer.',
                'fps': a.fps, 'poster_frame_index': len(frames)-1,
                'source': source, 'frames': frames}
    sequence_path.write_text(json.dumps(sequence, indent=2))
    review_path.write_text(json.dumps({'status': 'pending',
        'sequence_sha256': digest(sequence_path),
        'reason': 'Requires visual review of every actual frame before curator publication.'}, indent=2))
    print(sequence_path)


if __name__ == '__main__':
    main()
