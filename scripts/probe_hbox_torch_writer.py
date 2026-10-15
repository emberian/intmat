"""Resident three-array Torch counterpart of the opt-in Julia DoseThermalWriter.
Checks a Julia-exported multi-step fixture before timing one bounded grid.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import torch


def thermal_step(temperature, scratch, dose, ix, iy, iz, settings, source):
    dt, hx, hy, hz = settings[0], settings[1], settings[2], settings[3]
    x = ix*hx-settings[4]/2-source[0]
    y = iy*hy-settings[5]/2-source[1]
    z = iz*hz-settings[6]/2-source[2]
    heat = source[3]*torch.exp(-(x*x+y*y+z*z)/(2*settings[7]*settings[7]))
    lap = (torch.roll(temperature,1,0)+torch.roll(temperature,-1,0)-2*temperature)/(hx*hx)
    lap = lap+(torch.roll(temperature,1,1)+torch.roll(temperature,-1,1)-2*temperature)/(hy*hy)
    lap = lap+(torch.roll(temperature,1,2)+torch.roll(temperature,-1,2)-2*temperature)/(hz*hz)
    updated = torch.clamp(temperature+dt*(settings[8]*lap-settings[9]*temperature+heat),min=0)
    scratch.copy_(updated)
    dose.add_(dt*torch.clamp(updated-settings[10],min=0))


class ResidentWriter:
    def __init__(self, shape, fixture, operation):
        self.temperature = torch.zeros(shape,dtype=torch.float32,device='cuda')
        self.scratch = torch.zeros_like(self.temperature)
        self.dose = torch.zeros_like(self.temperature)
        self.axes = [(torch.arange(n,device='cuda',dtype=torch.float32)+.5).reshape(
            tuple(n if d==axis else 1 for d in range(3))) for axis,n in enumerate(shape)]
        h = [np.float32(e)/np.float32(n) for e,n in zip(fixture['extent'],shape)]
        dt = fixture['sources'][0]['dt']
        bound = .9/(fixture['cooling']+2*fixture['diffusion']*sum(1/float(v*v) for v in h))
        if dt > bound:
            raise ValueError(f'fixture dt {dt} exceeds thermal bound {bound} at {shape}')
        self.settings = torch.tensor([dt,*h,*fixture['extent'],fixture['width'],fixture['diffusion'],
            fixture['cooling'],fixture['threshold']],dtype=torch.float32,device='cuda')
        self.sources = torch.tensor([[*row['position'],row['power']] for row in fixture['sources']],
            dtype=torch.float32,device='cuda')
        if any(row['dt'] != dt for row in fixture['sources']):
            raise ValueError('this bounded probe requires fixed dt')
        self.operation = operation

    def advance(self, source):
        self.operation(self.temperature,self.scratch,self.dose,*self.axes,self.settings,source)
        self.temperature,self.scratch = self.scratch,self.temperature


def load_array(root, filename, shape):
    return np.fromfile(root/filename,dtype='<f4').reshape(shape,order='F')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('fixture')
    parser.add_argument('--n',type=int,default=384)
    parser.add_argument('--steps',type=int,default=25)
    parser.add_argument('--output',default='torch-writer-results.json')
    args=parser.parse_args()
    if not 16 <= args.n <= 512 or not 1 <= args.steps <= 100:
        raise ValueError('bounded grid 16..512 and steps 1..100 required')
    torch.set_num_threads(2)
    fixture_path=Path(args.fixture)
    fixture=json.loads(fixture_path.read_text())
    shape=tuple(fixture['shape'])
    operation=torch.compile(thermal_step,fullgraph=True)
    report=dict(torch=torch.__version__,hip=torch.version.hip,
        fixture_sha256=hashlib.sha256(fixture_path.read_bytes()).hexdigest(),fixture_source_hashes=fixture['sources_sha'],
        scope='opt-in three-array DoseThermalWriter; original explicit thermal integrator and analytic private material readout',
        frame_errors=[])
    started=time.perf_counter()
    writer=ResidentWriter(shape,fixture,operation)
    frames={frame['step']:frame for frame in fixture['frames']}
    for step,source in enumerate(writer.sources,1):
        writer.advance(source)
        if step in frames:
            row=dict(step=step)
            for field in ('temperature','dose'):
                actual=getattr(writer,field).cpu().numpy()
                expected=load_array(fixture_path.parent,frames[step][field],shape)
                error=float(np.max(np.abs(actual-expected)))
                row[field+'_max_abs_error']=error
                if not np.allclose(actual,expected,rtol=2e-5,atol=2e-6):
                    raise AssertionError(row)
            report['frame_errors'].append(row)
            print(json.dumps(row),flush=True)
    dose=writer.dose
    nutrient=fixture['initial_nutrient']*torch.exp(-fixture['conversion']*dose)
    deposit=fixture['initial_deposit']+fixture['initial_nutrient']*(-torch.expm1(-fixture['conversion']*dose))
    report['material_errors']={}
    for field,actual in [('nutrient',nutrient),('deposit',deposit)]:
        expected=load_array(fixture_path.parent,fixture['material_readout'][field],shape)
        values=actual.cpu().numpy()
        report['material_errors'][field+'_max_abs_error']=float(np.max(np.abs(values-expected)))
        assert np.allclose(values,expected,rtol=2e-5,atol=2e-6)
    report['fixture_compile_and_check_seconds']=time.perf_counter()-started
    del writer,dose,nutrient,deposit
    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats()
    writer=ResidentWriter((args.n,)*3,fixture,operation)
    compile_started=time.perf_counter()
    for source in writer.sources[:2]: writer.advance(source)
    torch.cuda.synchronize()
    report['benchmark_shape_compile_warm_seconds']=time.perf_counter()-compile_started
    writer.temperature.zero_(); writer.scratch.zero_(); writer.dose.zero_()
    torch.cuda.synchronize()
    started=time.perf_counter()
    for step in range(args.steps): writer.advance(writer.sources[step%len(writer.sources)])
    torch.cuda.synchronize()
    report['benchmark']=dict(grid=[args.n]*3,steps=args.steps,seconds=time.perf_counter()-started,
        dt=fixture['sources'][0]['dt'],source_schedule='fixture rows; repeats only if requested steps exceed fixture length',
        peak_allocated_vram=torch.cuda.max_memory_allocated(),peak_reserved_vram=torch.cuda.max_memory_reserved(),
        temperature_max=writer.temperature.max().item(),dose_max=writer.dose.max().item(),
        minimum_temperature=writer.temperature.min().item(),minimum_dose=writer.dose.min().item())
    report['benchmark']['seconds_per_step']=report['benchmark']['seconds']/args.steps
    assert report['benchmark']['minimum_temperature'] >= 0 and report['benchmark']['minimum_dose'] >= 0
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
    print('WRITER_PROBE_COMPLETE',flush=True)


if __name__=='__main__': main()
