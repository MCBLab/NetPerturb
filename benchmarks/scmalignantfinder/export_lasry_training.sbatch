#!/bin/bash
#SBATCH --partition=amd-3tb
#SBATCH --account=dmcoelho
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=300G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --job-name=export_lasry_train
#SBATCH --chdir=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder
#SBATCH --output=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/export_lasry_train_%j.out
#SBATCH --error=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/export_lasry_train_%j.err

set -euo pipefail
singularity exec --bind /mnt/beegfs --env R_MAX_VSIZE=300Gb /home/lgdqamorim/scratch/singularity_images/sctenifoldknk-v1.0.sif \
  Rscript /mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/export_lasry_training.R
