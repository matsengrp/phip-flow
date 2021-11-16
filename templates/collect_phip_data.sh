#!/bin/bash

set -euo pipefail

echo "Running phipflow"

phipflow load-from-counts-tsv \
    --sample_table !{sample_table} \
    --peptide_table !{peptide_table} \
    -c '*.counts' \
    -s '*.stats' \
    -o !{params.dataset_prefix}.phip

echo "Done"
