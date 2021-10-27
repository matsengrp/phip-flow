#!/bin/bash

set -euo pipefail

if [[ "${params.alignment_tool}" == "bowtie" ]]; then

    echo "Alignment tool is bowtie"

    ${params.fastq_stream_func} ${respective_replicate_path} | \
    bowtie ${params.align_args} --sam -x ${index}/peptide - > ${sample_id}.sam

else

    if [[ "${params.alignment_tool}" == "bowtie2" ]]; then

        echo "Alignment tool is bowtie2"

        ${params.fastq_stream_func} ${respective_replicate_path} | \
        bowtie2 ${params.align_args} -x ${index}/peptide - > ${sample_id}.sam

    else

        echo "${params.alignment_tool} is not recognized as bowtie or bowtie2 -- ERROR"

    fi

fi
