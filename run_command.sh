#!/bin/bash

nextflow run main.nf --obj ../data/vehicle_processed2.rds --column annotation --species mouse --n_cells 4000 --n_cores 4 --target ../alvos.txt --network genie3 -profile singularity -resume
