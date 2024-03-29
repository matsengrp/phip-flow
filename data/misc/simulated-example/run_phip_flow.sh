#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

/usr/bin/time nextflow  \
    -C ./nextflow.config \
    run ../main.nf \
    -with-report ./output/nextflow_report.html \
    -work-dir ./output/work/ \
    -resume
