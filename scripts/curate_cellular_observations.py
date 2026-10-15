"""Pack observed cellular geometry for the public living-material viewer.

Inputs are archived snapshots and explicitly labeled checkpoint continuations.
No simulation or generated/interpolated motion occurs here.
"""
import argparse
import hashlib
import json
from pathlib import Path


def pack(path):
    s = json.loads(path.read_text())
    agents = []
    for c in s['agents']:
        allocation = c['metabolism']['allocation']
        agents.append([c['id'], *c['position'][:2], c['radius'],
                       allocation[0] / max(sum(allocation[:2]), 1e-12),
                       c['energy'] / max(c['energy_capacity'], 1e-12),
                       c['metabolism']['adhesion'], c['generation'], c['age'],
                       c['hidden'][0], c['parent']])
    fields = {key: [round(x, 6) for x in s['fields'][key]['values']]
              for key in ('precursor_a', 'precursor_b')}
    m = s['metrics']
    return dict(step=s['step'], time=s['time'], cells=agents, bonds=m['cell_bonds'], fields=fields,
                metrics={k: m[k] for k in ('population', 'births', 'deaths', 'bonds', 'junction_transferred')},
                source_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('--destination', type=Path, default=Path('site/cellular'))
    a = p.parse_args(); a.destination.mkdir(parents=True, exist_ok=True)
    catalog = []
    for condition in ('broad', 'mosaic'):
        root = a.source / condition
        for kind in ('early', 'mature', 'history'):
            paths = sorted((root / ('continuation-' + kind)).glob('frame-*.json')) if kind != 'history' else sorted(root.glob('state-*.json'))
            assert paths, f'No snapshots for {condition}/{kind}'
            first = json.loads(paths[0].read_text())
            label = f'{condition.title()} · ' + {'early':'early contact network', 'mature':'later moving groups', 'history':'archived long trajectory'}[kind]
            note = ('Faithful continuation from an original saved checkpoint, with unchanged Julia code and parameters. These frames were reconstructed later; they were not recorded during the live soak.' if kind != 'history' else 'Actual archived states from the live soak. Large gaps between frames omit intervening movement and contact changes.')
            data = dict(label=label, note=note, grid=first['grid'], field_shape=first['fields']['precursor_a']['shape'],
                        columns=['id', 'x', 'y', 'radius', 'import_a_share', 'energy_fraction', 'adhesion', 'generation', 'age', 'recurrent_0', 'parent'],
                        frames=[pack(path) for path in paths])
            filename=f'{condition}-{kind}.json'
            (a.destination / filename).write_text(json.dumps(data, separators=(',', ':')))
            catalog.append(dict(label=label, file=filename, frames=len(paths), note=note))
    (a.destination / 'catalog.json').write_text(json.dumps(catalog, indent=2))
    print('Packed', sum(x['frames'] for x in catalog), 'frames; bytes', sum(p.stat().st_size for p in a.destination.glob('*.json')))


if __name__ == '__main__': main()
