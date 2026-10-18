#!/bin/bash
#SBATCH --partition=intel-128
#SBATCH --account=dmcoelho
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --job-name=scmalignantfinder_retrain_balanced
#SBATCH --chdir=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder
#SBATCH --output=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/retrain_bal_%j.out
#SBATCH --error=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/retrain_bal_%j.err

set -euo pipefail

module load softwares/python/3.10.5-gnu8
source /home/lgdqamorim/scratch/scmalignantfinder/venv/bin/activate
unset PYTHONPATH

# Training set is 64% Malignant / 36% Normal (Naldini+Beneyto combined), which
# biases the default LogisticRegression/RandomForest toward over-calling
# Malignant (confirmed: 92% sensitivity but only 51% specificity at the
# default 0.5 cutoff). init_core_model() in classifier.py doesn't expose
# class_weight, so patch it in our own venv copy to pass class_weight="balanced".
python - <<'PY'
p = "/home/lgdqamorim/scratch/scmalignantfinder/venv/lib/python3.10/site-packages/scMalignantFinder/classifier.py"
src = open(p).read()
old_lr = 'return LogisticRegression(n_jobs=n_thread)'
old_rf = 'return RandomForestClassifier(n_jobs=n_thread)'
new_lr = 'return LogisticRegression(n_jobs=n_thread, class_weight="balanced")'
new_rf = 'return RandomForestClassifier(n_jobs=n_thread, class_weight="balanced")'
changed = False
if old_lr in src:
    src = src.replace(old_lr, new_lr)
    changed = True
if old_rf in src:
    src = src.replace(old_rf, new_rf)
    changed = True
if changed:
    open(p, "w").write(src)
    print("patched init_core_model to use class_weight='balanced'")
else:
    print("no change needed (already patched or lines not found)")
PY

echo "=== training + predicting (LogisticRegression, class_weight=balanced) ==="
python predict_vangalen_retrained.py LogisticRegression logisticregression_balanced

echo "=== comparing ==="
python compare_predictions_generic.py vangalen_scmalignantfinder_predictions_retrained_logisticregression_balanced.csv

echo "=== done ==="
