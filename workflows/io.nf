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
process to_tall {

    label 'phippery'
    publishDir "$params.phip_data_dir/", mode: 'copy', overwrite: true

    input:
        file phip_data

    shell:
        """
        phippery to-tall-csv -o ${params.dataset_prefix}-tall.csv $phip_data 
        """
}

process to_wide {

    label 'phippery'
    publishDir "$params.phip_data_dir/", mode: 'copy', overwrite: true

    input:
        file phip_data

    shell:
        """
        mkdir $params.dataset_prefix
        phippery to-wide-csv -o $params.dataset_prefix $phip_data
        """
}
