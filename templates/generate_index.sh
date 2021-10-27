#!/bin/bash

set -euo pipefail

mkdir peptide_index

echo "Alignment tool is ${params.alignment_tool}"
echo "Running ${params.alignment_tool}-build"

${params.alignment_tool}-build \
    --threads ${task.cpus} \
    ${pep_fasta} \
    peptide_index/peptide
