#!/usr/bin/env python3
"""Build the dependency-free arena explorer beside an actual latest.json snapshot."""
import argparse
from pathlib import Path
import shutil

def build(destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(Path(__file__).resolve().parents[1] / 'web' / 'arena.html', destination / 'arena.html')
    return destination / 'arena.html'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination')
    print(build(parser.parse_args().destination))
