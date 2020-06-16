#!/bin/bash
set -e

nextflow \
    -C "nextflow.local.config" \
    run \
    main.nf \
    --params_file simulations/simulate_ones_single_ref.json \
    -with-report output/nextflow_report.html \
    -work-dir output/work/ \
    -ansi-log false \
