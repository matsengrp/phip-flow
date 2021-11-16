#!/bin/bash

set -euo pipefail

echo "Running phipflow"

phippery cpm -o !{phip_data} !{phip_data}
phippery query-expression "control_status=='library'" \
    -o lib.phip !{phip_data}
phippery fold-enrichment -dt "cpm" lib.phip !{phip_data}


# phipflow load-from-counts-tsv \
#     --sample_table !{sample_table} \
#     --peptide_table !{peptide_table} \
#     -c '*.counts' \
#     -s '*.stats' \
#     -o !{params.dataset_prefix}.phip
# 
echo "Done"
