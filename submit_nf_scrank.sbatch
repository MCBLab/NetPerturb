#!/bin/bash
#SBATCH --job-name=nf_scrank
#SBATCH --partition=intel-512
#SBATCH --account=dmcoelho
#SBATCH --qos=preempt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=480G
#SBATCH --exclude=r1i3n13
#SBATCH --time=2-00:00:00
#SBATCH --chdir=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank
#SBATCH --output=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/logs/nf_scrank_%j.out
#SBATCH --error=/mnt/beegfs/scratch/lgdqamorim/projeto_nextflow/nf_scrank/logs/nf_scrank_%j.err

set -euo pipefail

mkdir -p logs

# Remove o miniconda do PATH antes de iniciar o Nextflow.
# Sem isso, o Singularity herda o PATH do host e encontra o Python do miniconda
# (3.13) antes do Python interno do container, causando falha de codec/stdlib.
export PATH=$(echo "$PATH" | tr ':' '\n' | grep -v '/home/lgdqamorim/scratch/miniconda' | tr '\n' ':' | sed 's/:$//')
unset CONDA_DEFAULT_ENV CONDA_EXE CONDA_PYTHON_EXE CONDA_SHLVL CONDA_PREFIX CONDA_PROMPT_MODIFIER

# Informa ao Nextflow (executor local) os limites reais do nó alocado,
# evitando que ele tente lançar mais processos simultâneos do que o job suporta.
export NXF_OPTS="-Xms2g -Xmx8g"

nextflow run main.nf \
  --obj ../data/vehicle_processed2.rds \
  --column annotation \
  --species mouse \
  --n_cells 4000 \
  --n_cores 32 \
  --n_trees 50 \
  --target ../alvos.txt \
  --network genie3 \
  -profile singularity \
  -resume \
  -process.executor local \
  --max_cpus 32 \
  --max_memory 480.GB \
  --max_time 168.h
