"""Tiny extraction math gates; run with Blender --background --python this.py."""
import importlib.util
from pathlib import Path
import numpy as np
spec=importlib.util.spec_from_file_location('brine_surface',Path(__file__).with_name('build_brine_surface_explorer.py'));mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
# x-fast source positions and trilinear interpolation at native centers / boundaries.
shape=(4,3,2);field=np.fromfunction(lambda x,y,z:x+10*y+100*z,shape)
voxels=np.array([[0,0,0],[3,2,1],[1.5,.5,.25],[-1,0,0],[9,9,9]])
points=(voxels+.5)*2/np.array(shape)-1
np.testing.assert_allclose(mod.sample(field,points),[0,123,31.5,0,123],atol=1e-12)
values=np.array([[-1,3,0],[.123456,3,17],[1,3,100]])
encoded,lo,span,error=mod.quantize(values);restored=lo+encoded.astype(float)/65535*span
assert np.all(np.abs(values-restored)<np.array(span)/65535/2+1e-12)
assert encoded.dtype==np.dtype('<u2') and error[1]==0
print('brine surface math: x-fast native interpolation, edge clamp, quantization bounds passed',flush=True)
