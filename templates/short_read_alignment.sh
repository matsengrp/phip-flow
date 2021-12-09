#!/bin/bash

: '
This template aligns reads to the index after trimming the 
read to be the same length as tiles in the library.
For more on bowtie alignment, see the options
documented by bowtie.

We report only the best alignment found, and no more.

You may specify the number of allowed mismatches in the confi file
or hard-code it (or any other options) here.
'

set -euo pipefail

STREAM_FILE_CMD=!{params.fastq_stream_func}
FASTQ=!{respective_replicate_path}
INDEX=!{index}/peptide
ALIGN_OUT_FN=!{sample_id}.sam
READ_LENGTH=!{params.read_length}
PEPTIDE_LENGTH=!{params.peptide_tile_length}
CPUS=!{task.cpus}
MM=!{params.n_mismatches}

if [ ${PEPTIDE_LENGTH} -lt ${READ_LENGTH} ]; then
    let TRIM3=${READ_LENGTH}-${PEPTIDE_LENGTH}
else
    TRIM3=0
fi

: '
$STREAM_FILE_CMD $FASTQ | bowtie2 \
  -a \
  --trim3 $TRIM3 \
  --threads $CPUS \
  -x $INDEX - > $ALIGN_OUT_FN
'

$STREAM_FILE_CMD $FASTQ | bowtie \
  --trim3 $TRIM3 \
  --threads $CPUS \
  -n $MM \
  -l $PEPTIDE_LENGTH \
  -a \
  --tryhard \
  --nomaqround \
  --norc \
  --best \
  --sam \
  --quiet \
  -x $INDEX - > $ALIGN_OUT_FN









