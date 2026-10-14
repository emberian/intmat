"""Bounded JAX component probe; no births, exchanges, neighborhood search or arena step.

Float32 layout follows neural_controller!: [output,input,agent]. Gather supplies
six local field cues; remaining sensor cues are frozen fixture inputs. Physical
hidden/motion equations match neural_ecology.jl with externally supplied funding.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
import jax
import jax.numpy as jnp

H, S = 12, 32
DT, RESPONSE, SPEED = np.float32(.05), np.float32(1), np.float32(.3)
ALPHA = np.float32(-np.expm1(-float(DT) / float(RESPONSE)))
DIMS = np.array([32, 32, 1], dtype=np.float32)


def fixture(n):
    rng = np.random.default_rng(2718 + n)
    weights = rng.normal(0, .18, (H + 9, H + S, n)).astype('float32')
    marks = rng.normal(0, .01, weights.shape).astype('float32')
    hidden = rng.normal(0, .1, (H, n)).astype('float32')
    senses = rng.normal(0, .1, (S, n)).astype('float32')
    positions = (rng.random((3, n)) * DIMS[:, None]).astype('float32')
    positions[2] = 0
    fields = rng.random((6, 32, 32, 1)).astype('float32')
    drift = rng.normal(0, .01, (3, 32, 32, 1)).astype('float32')
    funding = rng.random(n).astype('float32')
    return weights, marks, hidden, senses, positions, fields, drift, funding


def component(weights, marks, hidden, senses, positions, fields, drift, funding):
    cells = jnp.floor(jnp.mod(positions, jnp.asarray(DIMS[:, None]))).astype(jnp.int32)
    local = fields[:, cells[0], cells[1], cells[2]]
    senses = senses.at[:6].set(local / (1 + jnp.abs(local)))
    inputs = jnp.concatenate((hidden, senses), axis=0)
    out = jnp.tanh(jnp.einsum('rcn,cn->rn', weights + marks, inputs,
                            precision=jax.lax.Precision.HIGHEST))
    hidden = hidden + (out[:H] - hidden) * ALPHA
    active = jnp.clip(DT * SPEED * out[H:H+3] * funding[None, :], -.45, .45)
    # Permeability is supplied in the sixth field cue, sampled at old/proposed cells.
    proposed = jnp.mod(positions + active, jnp.asarray(DIMS[:, None]))
    target = jnp.floor(proposed).astype(jnp.int32)
    permeability = jnp.minimum(local[5], fields[5, target[0], target[1], target[2]])
    flow = DT * drift[:, cells[0], cells[1], cells[2]]
    positions = jnp.mod(positions + jnp.clip(active * permeability[None, :] + flow, -.45, .45), jnp.asarray(DIMS[:, None]))
    positions = positions.at[2].set(0)
    return out, hidden, positions, senses


def reference(args):
    weights, marks, hidden, senses, positions, fields, drift, funding = args
    cells = np.floor(np.mod(positions, DIMS[:, None])).astype(np.int32)
    local = fields[:, cells[0], cells[1], cells[2]]
    senses = senses.copy(); senses[:6] = local / (1 + np.abs(local))
    inputs = np.concatenate((hidden, senses), axis=0)
    sums = np.zeros((H+9, len(funding)), dtype=np.float32)
    for col in range(H+S):
        sums += (weights[:, col, :] + marks[:, col, :]) * inputs[col]
    out = np.tanh(sums)
    hidden = hidden + (out[:H] - hidden) * ALPHA
    active = np.clip(DT*SPEED*out[H:H+3]*funding[None, :], -.45, .45)
    target = np.floor(np.mod(positions + active, DIMS[:, None])).astype(np.int32)
    permeability = np.minimum(local[5], fields[5, target[0], target[1], target[2]])
    positions = np.mod(positions + np.clip(active * permeability[None, :] + DT*drift[:, cells[0], cells[1], cells[2]], -.45, .45), DIMS[:, None])
    positions[2] = 0
    return out, hidden, positions, senses


def sync(result):
    return jax.block_until_ready(result)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--require-tpu', action='store_true')
    parser.add_argument('--output', default='tpu-controller-results.json')
    parser.add_argument('--sizes', default='64,1024,8192')
    args = parser.parse_args()
    report = {'scope': 'controller, local gather, physical update component only',
              'jax': jax.__version__, 'backend': jax.default_backend(),
              'devices': [str(d) for d in jax.devices()], 'measurements': []}
    print(json.dumps(report), flush=True)
    if args.require_tpu and report['backend'] != 'tpu':
        raise RuntimeError('TPU required; refusing CPU timings labeled as TPU')
    compiled = jax.jit(component)
    for n in map(int, args.sizes.split(',')):
        host = fixture(n)
        ref = reference(host)
        start = time.perf_counter(); device = jax.device_put(host); sync(device)
        transfer = time.perf_counter()-start
        start = time.perf_counter(); result = sync(compiled(*device))
        compile_first = time.perf_counter()-start
        errors = [float(np.max(np.abs(np.asarray(a)-b))) for a,b in zip(result,ref)]
        if max(errors) > 3e-5:
            raise AssertionError(f'Float32 CPU reference disagreement: {errors}')
        times = []
        for _ in range(20):
            start = time.perf_counter(); sync(compiled(*device)); times.append(time.perf_counter()-start)
        # Fused recurrent evolution excludes repeated host transfers and dispatch.
        def evolve(state):
            w,m,h,s,p,f,d,u = state
            def tick(carry, _):
                h,p,s=carry
                _,h,p,s=component(w,m,h,s,p,f,d,u)
                return (h,p,s), None
            return jax.lax.scan(tick,(h,p,s),None,length=100)[0]
        evolve_jit=jax.jit(evolve)
        start=time.perf_counter(); sync(evolve_jit(device)); recurrent_compile=time.perf_counter()-start
        start=time.perf_counter(); sync(evolve_jit(device)); recurrent=time.perf_counter()-start
        row={'agents': n, 'hidden': H, 'sensors': S, 'host_to_device_s': transfer,
             'compile_first_s': compile_first, 'warm_component_median_s': float(np.median(times)),
             'max_abs_errors_output_hidden_position_sensors': errors,
             'recurrent_100_compile_first_s': recurrent_compile, 'recurrent_100_warm_s': recurrent}
        report['measurements'].append(row); print(json.dumps(row),flush=True)
        Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print('PROBE_COMPLETE',flush=True)

if __name__ == '__main__':
    main()
