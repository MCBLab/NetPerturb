#!/bin/bash

nextflow run main.nf --obj ../data/GSE212217_pre_mutR_vs_NR.rds --column annotation --species human --n_cells 4000 --n_cores 32 --target ../alvos.txt --network genie3 -profile singularity -resume
