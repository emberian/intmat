"""Prepare a private, bounded Julia CUDA Kaggle kernel from the current source.

python3 scripts/prepare_kaggle_probe.py OUTPUT_DIRECTORY
Push the generated directory with Kaggle CLI and explicitly select NvidiaTeslaT4.
Only Project.toml, src/*.jl and the bounded benchmark are included.
"""
import base64
import io
import json
from pathlib import Path
import sys
import tarfile

root = Path(__file__).resolve().parents[1]
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
buffer = io.BytesIO()
with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
    for path in [root / "Project.toml", *sorted((root / "src").glob("*.jl")), root / "scripts/benchmark_accelerator.jl"]:
        archive.add(path, arcname=str(Path("intmat") / path.relative_to(root)))
payload = base64.b64encode(buffer.getvalue()).decode()
runner = '''import base64, json, os, pathlib, subprocess, tarfile, time, urllib.request
start=time.time()
output=pathlib.Path('/kaggle/working')
work=pathlib.Path('/tmp/intmat-probe')
work.mkdir(parents=True,exist_ok=True)
os.chdir(work)
subprocess.run(['nvidia-smi'],check=True)
(work/'source.tar.gz').write_bytes(base64.b64decode(PAYLOAD))
with tarfile.open(work/'source.tar.gz') as f: f.extractall(work, filter='data')
url='https://julialang-s3.julialang.org/bin/linux/x64/1.12/julia-1.12.4-linux-x86_64.tar.gz'
urllib.request.urlretrieve(url,work/'julia.tar.gz')
with tarfile.open(work/'julia.tar.gz') as f: f.extractall(work,filter='data')
julia=str(work/'julia-1.12.4/bin/julia')
env=dict(os.environ,JULIA_NUM_THREADS='2',JULIA_NUM_PRECOMPILE_TASKS='2',OPENBLAS_NUM_THREADS='1')
install='using Pkg; Pkg.activate("/tmp/intmat-probe/probe-env"); Pkg.develop(path="/tmp/intmat-probe/intmat"); Pkg.add("CUDA"); Pkg.precompile()'
subprocess.run([julia,'--startup-file=no','-e',install],env=env,check=True,timeout=1000)
print('INSTALL_SECONDS',time.time()-start,flush=True)
subprocess.run([julia,'--project=probe-env','-e','using CUDA; CUDA.versioninfo()'],env=env,check=True,timeout=90)
with open(output/'measurements.jsonl','w') as result:
 for backend in ['cpu','cuda']:
  for n in [32,64]:
   run=subprocess.run([julia,'--project=probe-env','intmat/scripts/benchmark_accelerator.jl',backend,str(n),'10'],env=env,check=True,timeout=150,text=True,capture_output=True)
   print(run.stdout,flush=True)
   print(run.stderr,flush=True)
   result.write(run.stdout); result.flush()
print('PROBE_COMPLETE',flush=True)
'''
(out / "run.py").write_text("PAYLOAD=" + repr(payload) + "\n" + runner)
(out / "kernel-metadata.json").write_text(json.dumps({
    "id": "emberian64/intmat-julia-cuda-bounded-probe", "title": "Intmat Julia CUDA bounded probe",
    "code_file": "run.py", "language": "python", "kernel_type": "script",
    "is_private": True, "enable_gpu": True, "enable_tpu": False,
    "enable_internet": True, "dataset_sources": [], "competition_sources": [],
    "kernel_sources": [], "model_sources": []}, indent=2) + "\n")
print(out)
