#!/bin/bash

set -euo pipefail

samtools stats ${sam_file} | \
    grep ^SN | 
    cut -f 2- | \
    sed '1p;7p;22p;25p;d' \
    > ${sample_id}.stats
