#!/bin/bash
#SBATCH --partition=intel-128
#SBATCH --account=dmcoelho
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --job-name=scmalignantfinder_retrain_xgb
#SBATCH --chdir=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder
#SBATCH --output=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/retrain_xgb_%j.out
#SBATCH --error=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/retrain_xgb_%j.err

set -euo pipefail

module load softwares/python/3.10.5-gnu8
source /home/lgdqamorim/scratch/scmalignantfinder/venv/bin/activate
# xgboost itself needs PYTHONPATH set so pip can run; keep it for the install.
pip install --use-deprecated=legacy-resolver xgboost
unset PYTHONPATH

echo "=== training + predicting (XGBoost, 3 cohorts) ==="
python predict_vangalen_retrained.py XGBoost xgboost_3cohorts

echo "=== comparing ==="
python compare_predictions_generic.py vangalen_scmalignantfinder_predictions_retrained_xgboost_3cohorts.csv

echo "=== done ==="
