#!/bin/bash
#SBATCH --partition=intel-128
#SBATCH --account=dmcoelho
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --job-name=vangalen_genie3
#SBATCH --chdir=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank
#SBATCH --output=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/slurm-%j.out

set -euo pipefail

# Van Galen run with the genie3 network method, same profile/flags as
# submit_vangalen_scrank.sbatch (which this was adapted from), to refresh the
# publication's Fig. 5 (fig2_vangalen.R, method=genie3): the previously
# published run's report predates the Data quality section (no identity_qc /
# network_qc / dq_flags tables), so this re-run is meant to replace it.

ROOT=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank
RUN=vangalen_genie3
RUNDIR="${ROOT}/nf_home/${RUN}"
OUTDIR="${ROOT}/results_genie3"
WORKDIR="${ROOT}/work_bench/${RUN}"
mkdir -p "$RUNDIR" "$WORKDIR"

export PATH=$(echo "$PATH" | tr ':' '\n' | grep -v '/home/lgdqamorim/scratch/miniconda' | tr '\n' ':' | sed 's/:$//')
unset CONDA_DEFAULT_ENV CONDA_EXE CONDA_PYTHON_EXE CONDA_SHLVL CONDA_PREFIX CONDA_PROMPT_MODIFIER

export NXF_OPTS="-Xms2g -Xmx8g"

cd "$RUNDIR"

nextflow run "${ROOT}/main.nf" \
  -profile singularity,test_vangalen \
  --outdir "${OUTDIR}" \
  --n_cores 32 \
  --network genie3 \
  --max_cpus 64 \
  --max_memory 64.GB \
  --max_time 12h \
  --min_cells 500 \
  -resume \
  --n_cells 2000 \
  --n_hvg 1000 \
  --target "${ROOT}/testdata/test_input_vangalen_expanded.txt" \
  -w "${WORKDIR}"
