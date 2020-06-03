#!/bin/bash
set -e

nextflow \
    -C "nextflow.local.config" \
    run \
    main.nf \
    --params_file simulations/simulate_ones_config.json \
    -with-report output/nextflow_report.html \
    -work-dir output/work/ \
    -ansi-log false \
    -with-dag dag.pdf \
    -resume
