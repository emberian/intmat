"""Package code-only private bounded TPU probe; token stays in Kaggle CLI.

python3 scripts/prepare_kaggle_tpu_probe.py OUTPUT_DIRECTORY
uvx --from kaggle kaggle kernels push -p OUTPUT_DIRECTORY --accelerator TpuV5E8 --timeout 600
"""
import json
from pathlib import Path
import sys

root=Path(__file__).resolve().parents[1]
out=Path(sys.argv[1]); out.mkdir(parents=True,exist_ok=True)
source=(root/'scripts/tpu_controller_probe.py').read_text()
(out/'probe.py').write_text(source)
# Require the TPU explicitly, without package upgrades or downloads.
(out/'run.py').write_text(source.replace('    main()\n', "    import sys\n    sys.argv.append('--require-tpu')\n    main()\n"))
(out/'kernel-metadata.json').write_text(json.dumps({
    'id':'emberian64/intmat-jax-tpu-controller-component-probe',
    'title':'Intmat JAX TPU controller component probe',
    'code_file':'run.py','language':'python','kernel_type':'script',
    'is_private':True,'enable_gpu':False,'enable_tpu':True,'enable_internet':False,
    'dataset_sources':[],'competition_sources':[],'kernel_sources':[],'model_sources':[]},indent=2)+'\n')
print(out)
