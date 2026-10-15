"""Bounded resident ROCm stencil/rFFT probe, independent of live arenas.
Run with the recovered hbox Torch interpreter and HSA_OVERRIDE_GFX_VERSION=10.3.0.
"""
import argparse
import json
import time
from pathlib import Path

import torch


def stencil(u):
    return u + 0.05 * (torch.roll(u, 1, 0) + torch.roll(u, -1, 0)
        + torch.roll(u, 1, 1) + torch.roll(u, -1, 1)
        + torch.roll(u, 1, 2) + torch.roll(u, -1, 2) - 6 * u)


def measure(fn, repetitions):
    fn()
    torch.cuda.synchronize()
    start = time.perf_counter()
    for _ in range(repetitions):
        fn()
    torch.cuda.synchronize()
    return (time.perf_counter() - start) / repetitions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sizes', default='64,128,256')
    parser.add_argument('--repetitions', type=int, default=5)
    parser.add_argument('--compile', action='store_true')
    parser.add_argument('--output', default='hbox-torch-results.json')
    args = parser.parse_args()
    sizes = [int(n) for n in args.sizes.split(',')]
    if not all(8 <= n <= 512 for n in sizes) or not 1 <= args.repetitions <= 20:
        raise ValueError('bounded probe requires sizes 8..512 and repetitions 1..20')
    torch.set_num_threads(2)
    report = dict(torch=torch.__version__, hip=torch.version.hip,
        device=torch.cuda.get_device_name(0), total_vram=torch.cuda.get_device_properties(0).total_memory,
        scope='periodic seven-point scalar diffusion and rFFT round trip; not Interior arena or complete writer',
        measurements=[])
    print(json.dumps(report), flush=True)
    reference = torch.linspace(-1, 1, 16**3).reshape(16, 16, 16)
    gpu = reference.cuda()
    parity = (stencil(gpu).cpu() - stencil(reference)).abs().max().item()
    assert parity < 1e-6, parity
    report['stencil_cpu_gpu_max_abs_error'] = parity
    compiled = None
    if args.compile:
        start = time.perf_counter()
        try:
            compiled = torch.compile(stencil, fullgraph=True)
            compiled_error = (compiled(gpu).cpu() - stencil(reference)).abs().max().item()
            assert compiled_error < 1e-6, compiled_error
            report['compiled_max_abs_error'] = compiled_error
            report['initial_compile_seconds'] = time.perf_counter() - start
        except Exception as error:
            compiled = None
            report['compile_error'] = f'{type(error).__name__}: {error}'
    for n in sizes:
        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats()
        u = torch.arange(n**3, dtype=torch.float32, device='cuda').reshape(n,n,n)
        u.remainder_(101).div_(101)
        host = torch.empty((n,n,n), dtype=torch.float32, pin_memory=True)
        host.copy_(u)
        destination = torch.empty_like(u)
        row = dict(n=n, scalar_bytes=u.numel()*u.element_size(),
            eager_stencil_seconds=measure(lambda: stencil(u), args.repetitions),
            host_to_device_seconds=measure(lambda: destination.copy_(host, non_blocking=True), args.repetitions),
            device_to_host_seconds=measure(lambda: host.copy_(u, non_blocking=True), args.repetitions))
        if compiled is not None:
            start = time.perf_counter()
            compiled(u)
            torch.cuda.synchronize()
            row['shape_compile_seconds'] = time.perf_counter() - start
            row['compiled_stencil_seconds'] = measure(lambda: compiled(u), args.repetitions)
        fft_start = time.perf_counter()
        spectrum = torch.fft.rfftn(u)
        restored = torch.fft.irfftn(spectrum, s=u.shape)
        torch.cuda.synchronize()
        row['fft_first_plan_and_roundtrip_seconds'] = time.perf_counter()-fft_start
        row['fft_roundtrip_max_abs_error'] = (restored-u).abs().max().item()
        assert row['fft_roundtrip_max_abs_error'] < 1e-5
        row['resident_rfft_roundtrip_seconds'] = measure(
            lambda: torch.fft.irfftn(torch.fft.rfftn(u), s=u.shape), args.repetitions)
        row['peak_allocated_vram'] = torch.cuda.max_memory_allocated()
        row['peak_reserved_vram'] = torch.cuda.max_memory_reserved()
        report['measurements'].append(row)
        Path(args.output).write_text(json.dumps(report, indent=2)+'\n')
        print(json.dumps(row), flush=True)
        del u, host, destination, spectrum, restored
    print('PROBE_COMPLETE', flush=True)


if __name__ == '__main__':
    main()
