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
process cpm_fold_enrichment {

    label 'phippery'

    input:
        file phip_data

    output:
        file phip_data

    shell:
        """
        phippery cpm -o !{phip_data} !{phip_data}
        phippery query-expression "control_status=='library'" \
            -o lib.phip !{phip_data}
        phippery fold-enrichment -dt "cpm" -o !{phip_data} lib.phip !{phip_data}
        """
}
