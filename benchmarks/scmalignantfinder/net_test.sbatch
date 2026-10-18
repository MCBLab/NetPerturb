#!/bin/bash
#SBATCH --partition=intel-128
#SBATCH --account=dmcoelho
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --job-name=net_test
#SBATCH --chdir=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder
#SBATCH --output=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/net_test_%j.out
#SBATCH --error=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/benchmarks/scmalignantfinder/net_test_%j.err
hostname
echo "--- pypi.org ---"
curl -sS -m 10 -o /dev/null -w "pypi http_code=%{http_code} time=%{time_total}s\n" https://pypi.org/simple/entrypoints/
echo "--- files.pythonhosted.org ---"
curl -sS -m 10 -o /dev/null -w "pythonhosted http_code=%{http_code} time=%{time_total}s\n" https://files.pythonhosted.org/
echo "--- zenodo (known to work) ---"
curl -sS -m 10 -o /dev/null -w "zenodo http_code=%{http_code} time=%{time_total}s\n" https://zenodo.org/
echo "done"
