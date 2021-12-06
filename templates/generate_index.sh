#!/bin/bash

set -euo pipefail

FASTA=!{pep_fasta}
CPUS=!{task.cpus}

mkdir peptide_index
bowtie2-build \
    --threads $CPUS \
    $FASTA \
    peptide_index/peptide
