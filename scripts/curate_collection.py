#!/usr/bin/env python3
"""Publish reviewed collection finals as 1600px WebP plus small metadata.

Run with a Python environment containing Pillow. Rendering is not performed.
Only final-status.json + matching art-reviewed.json can authorize an asset.
"""
import argparse
import hashlib
import json
import pathlib
import re
import sys

from build_pages import REPO, public


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked_path(root, value):
    path = pathlib.Path(value)
    if not path.is_absolute():
        # Renderer records either repository-relative or specimen-relative paths.
        path = REPO / path if (REPO / path).is_file() else root / path
    path = path.resolve()
    if root.resolve() not in path.parents:
        raise ValueError(f'Asset outside specimen: {value}')
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('batch', type=pathlib.Path)
    parser.add_argument('--ids', nargs='*', help='Optional subset of selected specimen IDs')
    args = parser.parse_args()
    batch = args.batch.resolve()
    if not re.fullmatch(r'collection-[0-9]{8}', batch.name):
        parser.error('batch must be a collection-YYYYMMDD directory')
    selected = set()
    for selection in batch.glob('collection-*-selected.json'):
        selected.update(read(selection))
    requested = set(args.ids) if args.ids else selected
    if requested - selected:
        parser.error('requested IDs must belong to the current selected manifests')
    from PIL import Image
    curated = []
    for slug in sorted(requested):
        if not re.fullmatch(r'[a-z0-9-]+', slug):
            raise ValueError(f'Invalid specimen ID: {slug}')
        origin = batch / slug
        final_path, review_path = origin / 'final-status.json', origin / 'art-reviewed.json'
        if not final_path.exists() or not review_path.exists():
            print(f'Skip {slug}: final render or final visual review pending')
            continue
        final, review = read(final_path), read(review_path)
        image = checked_path(origin, final.get('image', 'final.png'))
        source_hash = digest(image)
        if final.get('status') != 'complete' or final.get('stage') != 'final':
            raise ValueError(f'{slug}: not a completed final')
        if final.get('sha256') != source_hash or review.get('sha256') != source_hash or review.get('status') != 'approved':
            raise ValueError(f'{slug}: final hash/review mismatch')
        job = read(origin / 'render-job.json')
        if job.get('size') != 1600 or job.get('samples', 0) < 96:
            raise ValueError(f'{slug}: expected reviewed 1600px / 96-sample recipe')
        optical_path = origin / 'final.png.optical.json'
        if not optical_path.exists():
            raise ValueError(f'{slug}: optical metadata missing')
        with Image.open(image) as original:
            if original.width != 1600 or original.height != 1600:
                raise ValueError(f'{slug}: expected 1600px square final')
            target = REPO / 'site' / 'collection' / batch.name / slug
            target.mkdir(parents=True, exist_ok=True)
            original.convert('RGB').save(target / 'photograph.webp', 'WEBP', quality=95, method=6)
        run = public(read(origin / 'run.json'))
        optical_path = origin / 'final.png.optical.json'
        if not optical_path.exists():
            raise ValueError(f'{slug}: optical metadata missing')
        optical = public(read(optical_path))
        run['family'] = job.get('family', run.get('family', 'Material studies'))
        if not run.get('shape') and (origin / 'volume' / 'volume.json').exists():
            volume = read(origin / 'volume' / 'volume.json')
            run['shape'] = volume.get('shape', [])
        run['title'] = job.get('title', run.get('title', slug))
        run['description'] = job.get('description', run.get('description', ''))
        run['publication'] = dict(source_image_sha256=source_hash, webp_sha256=digest(target / 'photograph.webp'), quality=95, width=1600, original_retained=True, reviewed_by=review.get('reviewed_by', 'visual review'))
        for name, data in [('run.json', run), ('optical.json', optical)]:
            (target / name).write_text(json.dumps(data, indent=2) + '\n')
        recipe = public(read(origin / 'render-job.json'))
        (target / 'render-job.json').write_text(json.dumps(recipe, indent=2) + '\n')
        recipe_text = (origin / 'recipe.toml').read_text()
        if re.search(r'/Users/|/home/|file://|localhost:|127\.0\.0\.1:', recipe_text):
            raise ValueError(f'{slug}: recipe contains private references')
        (target / 'recipe.toml').write_text(recipe_text)
        for path in target.iterdir():
            if path.suffix in ('.json', '.toml') and re.search(r'/Users/|/home/|file://|localhost:|127\.0\.0\.1:', path.read_text()):
                raise ValueError(f'{slug}: private metadata reference')
        curated.append(dict(id=f'{batch.name}/{slug}', source_image_sha256=source_hash, files=sorted(p.name for p in target.iterdir() if p.is_file())))
        print(f'Curated {slug}: {(target / "photograph.webp").stat().st_size:,} bytes')
    # Extend publication history without removing previously published work.
    manifest_path = REPO / 'site' / 'collection-batches.json'
    existing = read(manifest_path) if manifest_path.exists() else []
    by_id = {item['id']: item for item in existing}
    by_id.update({item['id']: item for item in curated})
    if curated:
        manifest_path.write_text(json.dumps(list(by_id.values()), indent=2) + '\n')
    print(f'{len(curated)} reviewed finals curated; {len(requested) - len(curated)} pending')


if __name__ == '__main__':
    main()
