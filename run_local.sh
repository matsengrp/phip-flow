#!/bin/bash
set -e

nextflow \
    -C "nextflow.local.config" \
    run \
    main.nf \
    -params-file config.json \
    -with-report output/nextflow_report.html \
    -work-dir output/work/ \
    -resume
