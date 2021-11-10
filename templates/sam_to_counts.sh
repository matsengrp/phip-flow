#!/bin/bash

set -euo pipefail

# Convert SAM to BAM, and sort
samtools view -u -@ 28 ${sam_file} | \
    samtools sort -@ 28 - > ${sample_id}.bam

# Sort the BAM again
samtools sort -@ 28 ${sample_id}.bam -o ${sample_id}.sorted 

# Overwrite the first sorted BAM with the second
mv ${sample_id}.sorted ${sample_id}.bam

# Index the BAM
samtools index -b ${sample_id}.bam

# Count the number of reads per chromosome (excluding unmapped)
samtools idxstats ${sample_id}.bam | \
    cut -f 1,3 | \
    sed "/^*/d" > ${sample_id}.counts
