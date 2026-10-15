"""Native 3D field readback; accepts Python or Blender's -- arguments. No render."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from living_tissue_review import native, spectrum


def volume_review(root, output, times=(20.,40.,60.,80.)):
    root, output = Path(root), Path(output)
    if output.exists():
        raise FileExistsError(output)
    run = json.loads((root/'run.json').read_text()); dims = tuple(run['dimensions'])
    if len(dims)!=3 or min(dims)<=1:
        raise ValueError('native3D recording required')
    output.mkdir(); samples=[]
    for time in times:
        entry=min(run['history'],key=lambda x:abs(x['time']-time))
        directory=root/'frames'/f"{entry['step']:07d}"
        fields={name:native(directory,name,dims) for name in ('biomass','matrix','excitation','memory')}
        rho=fields['biomass']+fields['matrix']
        sample=dict(time=entry['time'],step=entry['step'],source_frame=str(directory),
                    quantile_levels=[.01,.1,.5,.9,.99],
                    quantiles={name:np.quantile(value,[.01,.1,.5,.9,.99]).tolist()
                               for name,value in (dict(fields,solid=rho)).items()},
                    low_solid_fractions={str(level):float((rho<level).mean()) for level in (.27,.45,.63)},
                    spectra={name:spectrum(value,run['extent']) for name,value in
                             [('excitation',fields['excitation']),('matrix',fields['matrix']),('solid',rho)]},
                    slices=[])
        for axis in range(3):
            index=(dims[axis]-1)//2
            target=output/f"{entry['step']:07d}"/'xyz'[axis];target.mkdir(parents=True)
            for name,field in fields.items():
                np.take(field,index,axis=axis).astype('<f4').ravel(order='F').tofile(target/(name+'.f32'))
            sample['slices'].append(dict(axis=axis,index=index,
                physical_coordinate=(index+.5)*run['extent'][axis]/dims[axis]-.5*run['extent'][axis],
                shape=[n for i,n in enumerate(dims) if i!=axis],directory=str(target.relative_to(output))))
        samples.append(sample)
    result=dict(format='living-tissue-native-volume-review-v1',source=str(root),
                dimensions=dims,extent=run['extent'],dt=run['dt'],capacity=run['parameters']['phase']['capacity'],
                source_sha256=run['source_sha256'],analysis_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in (Path(__file__),Path(__file__).with_name('living_tissue_review.py'))},samples=samples,
                convention='Exact quantiles and fluctuation spectra from complete native3D arrays. Slices are unresampled Float32 cell-center planes in Fortran order. Low solid fractions are B+M concentrations, not mechanical cavities or empty vacuum.')
    (output/'review.json').write_text(json.dumps(result,indent=2)+'\n')
    print('VOLUME REVIEW READY',output,flush=True)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root');parser.add_argument('output')
    args=parser.parse_args(sys.argv[sys.argv.index('--')+1:] if '--' in sys.argv else None)
    volume_review(args.root,args.output)
