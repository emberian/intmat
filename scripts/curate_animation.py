#!/usr/bin/env python3
"""Encode a reviewed sequence of rendered simulation frames for the gallery.

Input: sequence.json with id/title/description/fps/source and ordered frames
[{file,time,sha256}]. art-reviewed.json approves the SHA256 of sequence.json.
Times describe simulation samples; playback uses one video frame per sample.
The input directory and every original frame remain untouched.
"""
import argparse
import hashlib
import json
import math
import pathlib
import re
import shutil
import subprocess
import tempfile

from build_pages import REPO


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def public_metadata(value):
    if isinstance(value, dict):
        return {k: public_metadata(v) for k, v in value.items()}
    if isinstance(value, list):
        return [public_metadata(v) for v in value]
    if isinstance(value, str) and (value.startswith('/') or value.startswith('file://')):
        return pathlib.PurePosixPath(value).name
    return value


def validate_sequence(manifest):
    from PIL import Image
    manifest = manifest.resolve()
    origin = manifest.parent
    sequence = json.loads(manifest.read_text())
    review = json.loads((origin / 'art-reviewed.json').read_text())
    if review.get('status') != 'approved' or review.get('sequence_sha256') != digest(manifest):
        raise ValueError('Sequence must have a matching completed visual review')
    if not re.fullmatch(r'[a-z0-9]+(?:-[a-z0-9]+)*', sequence['id']):
        raise ValueError('Invalid animation id')
    fps = sequence['fps']
    if isinstance(fps, bool) or not isinstance(fps, (float, int)) or not math.isfinite(fps) or not 0 < fps <= 60:
        raise ValueError('fps must be finite and between zero and 60')
    frames = sequence['frames']
    if len(frames) < 2:
        raise ValueError('An animation needs at least two actual frames')
    poster_index = sequence.get('poster_frame_index', len(frames) // 2)
    if isinstance(poster_index, bool) or not isinstance(poster_index, int) or not 0 <= poster_index < len(frames):
        raise ValueError('Poster must select an actual reviewed frame index')
    if not sequence.get('source') or not sequence.get('title') or not sequence.get('description'):
        raise ValueError('Source, title and description are required')
    paths, dimensions, last_time = [], None, -math.inf
    for frame in frames:
        time = frame['time']
        if isinstance(time, bool) or not isinstance(time, (float, int)) or not math.isfinite(time) or time <= last_time:
            raise ValueError('Simulation times must be finite and strictly increasing')
        last_time = time
        path = (origin / frame['file']).resolve()
        if origin not in path.parents or path.suffix.lower() != '.png':
            raise ValueError('Frames must be PNG files inside the sequence directory')
        if digest(path) != frame['sha256']:
            raise ValueError(f'Frame changed after review: {path.name}')
        with Image.open(path) as image:
            image.verify()
        with Image.open(path) as image:
            if dimensions is None:
                dimensions = image.size
            elif image.size != dimensions:
                raise ValueError('Rendered frame dimensions must be constant')
        paths.append(path)
    return sequence, review, paths, dimensions


def curate(manifest, destination, *, max_width=1280, budget_bytes=12_000_000):
    from PIL import Image
    if max_width < 2 or budget_bytes <= 0:
        raise ValueError('Positive artifact budget and width >= 2 required')
    sequence, review, paths, dimensions = validate_sequence(manifest)
    manifest_hash = review['sequence_sha256']
    target = destination.resolve() / sequence['id']
    if target.exists():
        raise ValueError('Publication exists; select a new version id')
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.animation-', dir=target.parent) as temporary:
        work = pathlib.Path(temporary)
        frames = work / 'frames'
        frames.mkdir()
        for index, path in enumerate(paths):
            # Symlinks avoid copying potentially gigabytes of reviewed frames.
            (frames / f'{index:06d}.png').symlink_to(path)
        video = work / 'film.mp4'
        subprocess.run([
            'ffmpeg', '-hide_banner', '-loglevel', 'error', '-nostdin',
            '-framerate', str(sequence['fps']), '-i', str(frames / '%06d.png'),
            '-frames:v', str(len(paths)), '-an', '-c:v', 'libx264', '-preset', 'slow',
            '-crf', '18', '-threads', '2', '-pix_fmt', 'yuv420p',
            '-vf', f"scale='min(iw,{max_width})':-2:flags=lanczos,pad=ceil(iw/2)*2:ceil(ih/2)*2",
            '-movflags', '+faststart', str(video),
        ], check=True)
        probe = json.loads(subprocess.check_output([
            'ffprobe', '-v', 'error', '-count_frames', '-select_streams', 'v:0',
            '-show_entries', 'stream=width,height,nb_read_frames,duration', '-of', 'json', str(video)
        ]))['streams'][0]
        if int(probe['nb_read_frames']) != len(paths):
            raise ValueError('Encoded frame count differs from reviewed sequence')
        expected_duration = len(paths) / sequence['fps']
        if abs(float(probe['duration']) - expected_duration) > .02:
            raise ValueError('Encoded duration differs from frame-time mapping')
        if digest(manifest) != manifest_hash or any(digest(path) != frame['sha256']
                for path, frame in zip(paths, sequence['frames'])):
            raise ValueError('Reviewed source changed during encoding')
        poster_index = sequence.get('poster_frame_index', len(paths) // 2)
        with Image.open(paths[poster_index]) as image:
            image.convert('RGB').save(work / 'poster.webp', 'WEBP', quality=93, method=6)
        public = {
            'schema': 'interior-animation', 'id': sequence['id'],
            'title': sequence['title'], 'description': sequence['description'],
            'source': public_metadata(sequence['source']),
            'playback': {'fps': sequence['fps'], 'duration_seconds': expected_duration,
                         'width': probe['width'], 'height': probe['height'],
                         'rendered_width': dimensions[0], 'rendered_height': dimensions[1],
                         'interpolation': 'none; one video frame per recorded simulation sample'},
            'frames': [{'video_time': i / sequence['fps'], 'simulation_time': frame['time'],
                        'source_image_sha256': frame['sha256']} for i, frame in enumerate(sequence['frames'])],
            'review': {'sequence_sha256': manifest_hash, 'reviewed_by': review.get('reviewed_by', 'visual review')},
            'poster_frame_index': poster_index,
            'video_sha256': digest(video), 'poster_sha256': digest(work / 'poster.webp'),
        }
        (work / 'film.json').write_text(json.dumps(public, indent=2) + '\n')
        files = ('film.mp4', 'film.json', 'poster.webp')
        size = sum((work / name).stat().st_size for name in files)
        if size > budget_bytes:
            raise ValueError(f'Animation exceeds publication budget: {size} > {budget_bytes}')
        # Atomic directory publication after encoding, decode-count validation and budgeting.
        staged = work / 'published'
        staged.mkdir()
        for name in files:
            shutil.move(work / name, staged / name)
        staged.rename(target)
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=pathlib.Path)
    parser.add_argument('--destination', type=pathlib.Path, default=REPO / 'site' / 'animations')
    parser.add_argument('--max-width', type=int, default=1280)
    parser.add_argument('--budget-bytes', type=int, default=12_000_000)
    args = parser.parse_args()
    print(curate(args.manifest, args.destination, max_width=args.max_width, budget_bytes=args.budget_bytes))


if __name__ == '__main__':
    main()
