#!/usr/bin/env nextflow
/*
Compute fold enrichment workflow

Author: Jared G. Galloway
*/

// Using DSL-2
nextflow.enable.dsl=2

// TODO eventually, we'll just make a generic workflow for doing
// all the types of normalization. diff sel, model fitting (p-val), etc.
// all that layer the analysis.

// compute fold enrichment on a single .phip file
process fold_enrichment_process {

    label 'phippery'

    input:
        file phip_data

    output:
        file phip_data

    shell:
    template fold_enrichment.sh
}

workflow FOLD_ENR {
    take: phip_data
    main:
        fold_enrichment_process(phip_data)
    emit:
        fold_enrichment_process.out
}
