#!/bin/bash
#SBATCH --partition=intel-128
#SBATCH --account=dmcoelho
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --job-name=vangalen_sctknk_isolated
#SBATCH --chdir=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank
#SBATCH --output=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/slurm-%j.out

set -euo pipefail

# Same vangalen_sctknk run as submit_vangalen_sctknk.sbatch, but with
# --network forced to null (conf/no_network.config, which overrides the
# test_vangalen profile's network="hdwgcna" default) so only the --sctknk
# track runs: no HDWGCNA, no RANK_SCORE/MERGE/NETWORK_QC/REPORT.
#
# Same RUN name and -w work dir as the combined run on purpose: DOWNSAMPLE
# and the whole SCTENIFOLDKNK_BUILD/NULL/KO/MERGE_SCTENIFOLDKNK track never
# took `network` as an input, so their task hashes are identical either
# way -- -resume should pull all of that straight from cache instead of
# recomputing the 1045 KO tasks. Only OUTDIR differs, so the (nonexistent,
# with network=null) report/rank-score outputs don't collide with the
# combined run's published ones.

ROOT=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank
RUN=vangalen_sctknk
RUNDIR="${ROOT}/nf_home/${RUN}"
OUTDIR="${ROOT}/results_sctknk_isolated"
WORKDIR="${ROOT}/work_bench/${RUN}"
mkdir -p "$RUNDIR" "$WORKDIR" "$OUTDIR"

export PATH=$(echo "$PATH" | tr ':' '\n' | grep -v '/home/lgdqamorim/scratch/miniconda' | tr '\n' ':' | sed 's/:$//')
unset CONDA_DEFAULT_ENV CONDA_EXE CONDA_PYTHON_EXE CONDA_SHLVL CONDA_PREFIX CONDA_PROMPT_MODIFIER

export NXF_OPTS="-Xms2g -Xmx8g"

cd "$RUNDIR"

nextflow run "${ROOT}/main.nf" \
  -profile singularity,test_vangalen \
  -c "${ROOT}/conf/no_network.config" \
  --outdir "${OUTDIR}" \
  --n_cores 32 \
  --sctknk true \
  --max_cpus 64 \
  --max_memory 64.GB \
  --max_time 12h \
  --min_cells 500 \
  -resume \
  --n_cells 2000 \
  --n_hvg 1000 \
  --target "${ROOT}/testdata/test_input_vangalen_expanded.txt" \
  -w "${WORKDIR}"
