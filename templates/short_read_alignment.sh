#!/bin/bash

set -euo pipefail

STREAM_FILE_CMD=!{params.fastq_stream_func}
FASTQ=!{respective_replicate_path}
INDEX=!{index}/peptide
ALIGN_OUT_FN=!{sample_id}.sam

$STREAM_FILE_CMD $FASTQ | bowtie \
  --trim3 8 \
  --threads 4 \
  -n 2 \
  -l 117 \
  --tryhard \
  --nomaqround \
  --norc \
  --best \
  --sam \
  --quiet \
  -x $INDEX - > $ALIGN_OUT_FN
