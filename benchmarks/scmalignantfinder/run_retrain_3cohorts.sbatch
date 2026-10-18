#!/bin/bash
#SBATCH --partition=intel-128
#SBATCH --account=dmcoelho
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --job-name=scmalignantfinder_retrain_3cohorts
#SBATCH --chdir=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder
#SBATCH --output=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/retrain_3c_%j.out
#SBATCH --error=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/retrain_3c_%j.err

set -euo pipefail

module load softwares/python/3.10.5-gnu8
source /home/lgdqamorim/scratch/scmalignantfinder/venv/bin/activate
unset PYTHONPATH

echo "=== building training h5ad (Naldini + Beneyto + Lasry) ==="
python build_train_h5ad.py

echo "=== training + predicting (LogisticRegression) ==="
python predict_vangalen_retrained.py LogisticRegression logisticregression_3cohorts

echo "=== comparing ==="
python compare_predictions_generic.py vangalen_scmalignantfinder_predictions_retrained_logisticregression_3cohorts.csv

echo "=== done ==="
