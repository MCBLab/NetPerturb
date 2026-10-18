#!/bin/bash
#SBATCH --partition=amd-3tb
#SBATCH --account=dmcoelho
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --job-name=scmalignantfinder_annotate
#SBATCH --chdir=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder
#SBATCH --output=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/annotate_%j.out
#SBATCH --error=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/annotate_%j.err

set -euo pipefail

echo "=== exporting all cells from the filtered pipeline objects ==="
singularity exec --bind /mnt/beegfs --env R_MAX_VSIZE=200Gb /home/lgdqamorim/scratch/singularity_images/sctenifoldknk-v1.0.sif \
  Rscript /mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/export_for_annotation.R

module load softwares/python/3.10.5-gnu8
source /home/lgdqamorim/scratch/scmalignantfinder/venv/bin/activate
unset PYTHONPATH

for study in Naldini_2023 BeneytoCalabuig_2023 Lasry_2022; do
  echo "=== predicting: $study ==="
  python predict_annotation.py "$study"
done

module unload softwares/python/3.10.5-gnu8 2>/dev/null || true
echo "=== merging predictions back into the filtered .rds objects ==="
singularity exec --bind /mnt/beegfs --env R_MAX_VSIZE=200Gb /home/lgdqamorim/scratch/singularity_images/sctenifoldknk-v1.0.sif \
  Rscript /mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/merge_predictions_into_rds.R

echo "=== done ==="
