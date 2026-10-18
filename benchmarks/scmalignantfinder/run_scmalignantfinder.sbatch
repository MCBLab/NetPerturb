#!/bin/bash
#SBATCH --partition=intel-128
#SBATCH --account=dmcoelho
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --job-name=scmalignantfinder_vangalen
#SBATCH --chdir=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder
#SBATCH --output=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/run_%j.out
#SBATCH --error=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/run_%j.err

set -euo pipefail

module load softwares/python/3.10.5-gnu8

VENV=/home/lgdqamorim/scratch/scmalignantfinder/venv
rm -rf "$VENV"
python3 -m venv "$VENV"
source "$VENV/bin/activate"
# This venv has no self-contained pip (ensurepip didn't bootstrap one); it
# only works via the module's PYTHONPATH, so keep it for every `pip` call.
# But the module's PYTHONPATH also sits ahead of the venv's own site-packages
# in sys.path, so anything installed into the venv (pandas 2.1.4, say) is
# silently shadowed by the module's older copy (pandas 1.4.3) at import time.
# Fix: unset PYTHONPATH only around running our own scripts, not around pip.

export PIP_DEFAULT_TIMEOUT=120
python -m pip install --upgrade pip -q
# The new resolver backtracks badly over scMalignantFinder's pinned heavy deps
# (pyscenic, squidpy, dask) and can run for hours; every dep here is already
# an exact pin, so the legacy resolver (first-match, no cross-graph search)
# is both faster and sufficient. Install ONLY scMalignantFinder: adding
# unpinned anndata/scipy/pandas here outranks its own exact pins under the
# legacy resolver and leaves an inconsistent environment (seen in job 2206982:
# pandas 1.4.3 instead of the required 2.1.4, breaking dask/spatialdata).
# --ignore-installed forces every dependency into the venv itself: with
# PYTHONPATH pointing at the module's global site-packages, pip otherwise
# treats small packages already there (h5py, requests, jinja2, six, ...) as
# "already satisfied" and never copies them into the venv -- fine while
# PYTHONPATH stays set, but they vanish once we unset it to run our scripts
# (needed so the venv's own pandas/numpy/anndata aren't shadowed instead).
# setuptools>=81 dropped the pkg_resources shim that pyscenic/ctxcore/loompy
# (old-style packages) still `import pkg_resources` directly; pin it back.
pip install --use-deprecated=legacy-resolver --ignore-installed --retries 10 "setuptools<81" scMalignantFinder

echo "=== packages installed, running prediction ==="
# scMalignantFinder/__init__.py unconditionally does `from . import spatial`,
# which pulls in squidpy -> spatialdata -> ome-zarr, and that chain needs
# zarr>=3 (FsspecStore) while something else in the tree pins zarr 2.x --
# an upstream conflict in a module (spatial region ID) we don't use at all;
# classifier.py itself only needs numpy/pandas/scanpy/scipy/sklearn. Make
# that import optional in our own venv copy rather than fight zarr/squidpy
# version hell for a feature this run never calls.
python - <<'PY'
import re
p = "/home/lgdqamorim/scratch/scmalignantfinder/venv/lib/python3.10/site-packages/scMalignantFinder/__init__.py"
src = open(p).read()
patched = src.replace(
    "from . import spatial",
    "try:\n    from . import spatial\nexcept ImportError:\n    pass",
)
if patched != src:
    open(p, "w").write(patched)
    print("patched __init__.py to make `spatial` import optional")
else:
    print("no change needed (already patched or line not found)")
PY

unset PYTHONPATH
python predict_vangalen.py

echo "=== running comparison ==="
python compare_predictions.py

echo "=== done ==="
