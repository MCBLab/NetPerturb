#!/bin/bash
#SBATCH --partition=intel-128
#SBATCH --account=dmcoelho
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --job-name=scmalignantfinder_retrain_rf
#SBATCH --chdir=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder
#SBATCH --output=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/retrain_rf_%j.out
#SBATCH --error=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/retrain_rf_%j.err

set -euo pipefail

module load softwares/python/3.10.5-gnu8
source /home/lgdqamorim/scratch/scmalignantfinder/venv/bin/activate
unset PYTHONPATH

# Reuses combine_training_aml.h5ad (Naldini+Beneyto+Lasry, 295k cells),
# already built for the LogisticRegression run -- no need to rebuild it.
echo "=== training + predicting (RandomForest, 3 cohorts) ==="
python predict_vangalen_retrained.py RandomForest randomforest_3cohorts

echo "=== comparing ==="
python compare_predictions_generic.py vangalen_scmalignantfinder_predictions_retrained_randomforest_3cohorts.csv

echo "=== done ==="
