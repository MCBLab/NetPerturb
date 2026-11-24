#!/bin/bash
#SBATCH --partition=intel-128
#SBATCH --account=dmcoelho
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=00:20:00
#SBATCH --job-name=fig2_vangalen
#SBATCH --chdir=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/NetPerturb/publication/analysis
#SBATCH --output=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/NetPerturb/publication/analysis/fig2_vangalen_job_%j.out

set -euo pipefail

METHOD="${1:-genie3}"

singularity exec --bind /mnt/beegfs \
  /home/lgdqamorim/scratch/singularity_images/juliaapolonio-scrank-latest.img \
  Rscript fig2_vangalen.R "$METHOD"
