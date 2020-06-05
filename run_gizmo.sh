#!/bin/bash
set -e
source /app/Lmod/lmod/lmod/init/bash
module use /app/easybuild/modules/all

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

nextflow \
    -C "/home/jgallowa/.nextflow/config" \
    run \
    main_simulation.nf \
    --params_file simulations/simulate_ones_config.json \
    -with-report output/nextflow_report.html \
    -work-dir /fh/scratch/delete30/matsen_e/jgallowa/temp/work/ \
    -ansi-log false \
    -resume
