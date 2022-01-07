#!/usr/bin/env nextflow
/*
Compute fold enrichment workflow

Author: Jared G. Galloway
*/

// Using DSL-2
nextflow.enable.dsl=2

/*
AUTOMATICALLY COMPUTED
----------------------
(NO REQUIRED ANNOTATIONS)
*/

process counts_per_million {
    input: path phip_data
    output: path "cpm.phip"
    shell:
    """
    phippery cpm -o cpm.phip ${phip_data} 
    """
}

process size_factors {
    input: path phip_data
    output: path "sf.phip"
    shell:
    """
    phippery size-factors -o sf.phip ${phip_data} 
    """
}

// TODO
//process rank_counts {
//    input:
//    path phip_data
//    output:
//    path "rank.phip"
//    shell:
//    """
//    phippery rank -o rank.phip !{phip_data} 
//    """
//}


/*
OPTIONALLY RUN STATISTICS
-------------------------
(ANNOTATIONS REQUIRED & FLAG)
*/

process cpm_fold_enrichment {
    input: path phip_data
    output: path "fold_enr.phip"
    when: params.run_cpm_enr_workflow
    shell:
    """
    phippery query-expression "control_status=='library'" \
        -o lib.phip ${phip_data}
    phippery fold-enrichment -dt "cpm" -o fold_enr.phip lib.phip ${phip_data}
    """
}

// compute fold enrichment on a single .phip file
//process differential_selection {
//    // label "phip_dev"
//    // conda = '/home/jared/miniconda3/envs/phip_seq/'
//    input:
//    path phip_data
//    output:
//    path "diff_sel.phip"
//    when: params.run_diff_sel_workflow
//    shell:
//    """
//    cp ${phip_data} diff_sel.phip
//    """
//}


/*
a generic process using the xarray
merge infrastructure
*/

process merge_binary_datasets {    
    input:
    path all_phip_datasets
    output:
    path "merged.phip"
    shell:
    """
    phippery merge -o merged.phip '*.phip'
    """
}


workflow STATS {

    take: dataset
    main:

    // we automatically compute some stats
    // which are independent of any annotations
    dataset | \
        (counts_per_million & size_factors) | \
        mix | set { auto_stats_ch }

    // run some optional statistics which
    // depend on certain annotations
    cpm_fold_enrichment(counts_per_million.out) | set { cpm_fold_enr }

    // TODO diff sel
    // TODO neg binom

    // collect all the datasets statistics and merge
    auto_stats_ch.concat(
        cpm_fold_enr
    ) | collect | merge_binary_datasets

    emit:
    merge_binary_datasets.out
} 
