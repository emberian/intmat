#!/usr/bin/env python3
"""Publish reviewed collection finals as full-resolution WebP plus small metadata.

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
    if not re.fullmatch(r'(collection|refinement)-[0-9]{8}', batch.name):
        parser.error('batch must be collection-YYYYMMDD or refinement-YYYYMMDD')
    selected = {}
    explicit = batch / 'refinement-publication-selected.json'
    selections = [explicit] if batch.name.startswith('refinement-') and explicit.exists() else list(batch.glob('collection-*-selected.json')) if batch.name.startswith('collection-') else []
    for selection in selections:
        for item in read(selection):
            record = {'id': item} if isinstance(item, str) else item
            selected[record['id']] = record
    requested = set(args.ids) if args.ids else set(selected)
    if requested - set(selected):
        parser.error('requested IDs must belong to the current selected manifests')
    from PIL import Image
    curated = []
    for slug in sorted(requested):
        if not re.fullmatch(r'[a-z0-9-]+', slug):
            raise ValueError(f'Invalid specimen ID: {slug}')
        record = selected[slug]
        origin = pathlib.Path(record.get('directory', slug))
        if not origin.is_absolute():origin = REPO / origin if (REPO / origin).is_dir() else batch / origin
        origin = origin.resolve()
        if batch not in origin.parents:raise ValueError(f'{slug}: directory outside batch')
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
        size, samples = job.get('size', 0), job.get('samples', 0)
        minimum_samples = record.get('minimum_samples', 96)
        if type(minimum_samples) is not int or minimum_samples < 1:
            raise ValueError(f'{slug}: minimum_samples must be a positive integer')
        sampling_review = review.get('sampling_review')
        if minimum_samples < 96 and (not isinstance(sampling_review, str) or not sampling_review.strip()):
            raise ValueError(f'{slug}: a lower sample budget needs a recorded visual sampling review')
        if size < 1600 or samples < minimum_samples:
            raise ValueError(f'{slug}: expected at least 1600px / {minimum_samples}-sample recipe')
        for key in ('size','samples'):
            if key in record and record[key] != job[key]:raise ValueError(f'{slug}: selection and job {key} differ')
        run = public(read(origin / 'run.json'))
        volume_path = origin / job.get('volume', 'volume') / 'volume.json'
        volume = read(volume_path) if volume_path.exists() else {}
        native = record.get('native_shape', run.get('native_shape', []))
        if batch.name.startswith('refinement-'):
            if not native or native != run.get('shape', volume.get('shape')):
                raise ValueError(f'{slug}: declared native shape does not match recorded source')
            run['native_shape'] = native
            run['comparison_notes'] = record.get('comparison_notes', run.get('comparison_notes', 'Native solver grid; photograph pixel dimensions are listed separately.'))
        optical_path = checked_path(origin, record.get('optical', final.get('optical', 'final.png.optical.json')))
        if not optical_path.exists():
            raise ValueError(f'{slug}: optical metadata missing')
        optical = public(read(optical_path))
        settings = optical.get('settings', {})
        if settings.get('samples') != samples or settings.get('size') != size:
            raise ValueError(f'{slug}: optical render settings do not match selected job')
        with Image.open(image) as original:
            if original.width != size or original.height != size:
                raise ValueError(f'{slug}: final dimensions do not match {size}px job')
            target = REPO / 'site' / 'collection' / batch.name / slug
            target.mkdir(parents=True, exist_ok=True)
            original.convert('RGB').save(target / 'photograph.webp', 'WEBP', quality=95, method=6)
        run['family'] = job.get('family', run.get('family', 'Material studies'))
        if not run.get('shape') and (origin / 'volume' / 'volume.json').exists():
            volume = read(origin / 'volume' / 'volume.json')
            run['shape'] = volume.get('shape', [])
        run['title'] = job.get('title', run.get('title', slug))
        run['description'] = re.sub(r'Native(?=\d)', 'Native ', job.get('description', run.get('description', ''))).replace('endtime', 'end time')
        if run.get('comparison_notes'):run['comparison_notes'] = re.sub(r'Native(?=\d)', 'Native ', run['comparison_notes'])
        run['publication'] = dict(source_image_sha256=source_hash, webp_sha256=digest(target / 'photograph.webp'), quality=95, width=size, samples=samples, native_shape=native, original_retained=True, reviewed_by=review.get('reviewed_by', 'visual review'))
        if minimum_samples < 96:
            run['publication']['sampling_review'] = sampling_review.strip()
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
