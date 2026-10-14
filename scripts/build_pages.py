#!/usr/bin/env python3
"""Stage the curated, tracked gallery; --refresh updates it from local renders."""
import argparse
import json
import pathlib
import re
import shutil
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
SOURCE = REPO / 'site' / 'collection'
MANIFEST = REPO / 'site' / 'manifest.json'
GITHUB = 'https://github.com/emberian/intmat/blob/dev/'

def public(value):
    if isinstance(value, dict):
        return {k: public(v) for k, v in value.items()}
    if isinstance(value, list):
        return [public(v) for v in value]
    if isinstance(value, str) and (value.startswith(('/Users/', '/home/', '/tmp/')) or value.startswith('file://')):
        return pathlib.Path(value).name
    return value

def refresh():
    selection = json.loads(MANIFEST.read_text())
    # Only replace files in the dedicated generated collection, never outputs.
    for item in selection:
        origin = REPO / 'outputs' / item['source']
        target = SOURCE / item['id']
        target.mkdir(parents=True, exist_ok=True)
        for old in target.iterdir():
            if old.is_file():
                old.unlink()
        for name in item['files']:
            source = origin / name
            if not source.is_file():
                continue
            dest = target / name
            if source.suffix == '.json':
                dest.write_text(json.dumps(public(json.loads(source.read_text())), indent=2) + '\n')
            elif source.suffix == '.html':
                text = source.read_text()
                def clean(match):
                    return match[1] + json.dumps(public(json.loads(match[2])), separators=(',', ':')).replace('</', '<\\/') + match[3]
                text = re.sub(r'(<script id="recording" type="application/json">)(.*?)(</script>)', clean, text, flags=re.S)
                dest.write_text(text)
            else:
                shutil.copyfile(source, dest)
        # Optical-only recordings use volume metadata as their gallery manifest.
        volume = origin / 'volume.json'
        if not volume.exists():
            volume = origin / 'volume' / 'volume.json'
        if not (target / 'run.json').exists() and volume.exists():
            (target / 'volume.json').write_text(json.dumps(public(json.loads(volume.read_text())), indent=2) + '\n')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh', action='store_true')
    parser.add_argument('--destination', default='_site')
    args = parser.parse_args()
    if args.refresh:
        refresh()
    destination = pathlib.Path(args.destination).resolve()
    if destination == SOURCE or SOURCE in destination.parents or destination == REPO:
        parser.error('destination must be separate from the source collection and repository root')
    destination.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SOURCE, destination, dirs_exist_ok=True)
    subprocess.run([sys.executable, str(REPO / 'scripts/build_gallery.py'), str(destination)], check=True)
    index = destination / 'index.html'
    text = index.read_text()
    # Documentation lives in the repository, independent of the Pages base path.
    text = text.replace('../docs/', GITHUB + 'docs/').replace('local collection', 'recorded material collection')
    text = text.replace('</header>', '<a href="https://github.com/emberian/intmat">Source &amp; simulation code ↗</a></header>')
    index.write_text(text)
    for file in destination.rglob('*'):
        if file.is_file() and file.suffix in ('.html', '.json', '.toml'):
            if re.search(r'/Users/|/home/|file://|localhost:|127\.0\.0\.1:', file.read_text()):
                raise ValueError(f'Private/local reference in {file.relative_to(destination)}')
    size = sum(p.stat().st_size for p in destination.rglob('*') if p.is_file())
    if size > 30_000_000:
        raise ValueError(f'Gallery exceeds 30 MB budget: {size}')
    print(f'Staged {size / 1_000_000:.2f} MB at {destination}')

if __name__ == '__main__':
    main()
