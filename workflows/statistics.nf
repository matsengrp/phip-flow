#!/usr/bin/env nextflow
/*
Compute fold enrichment workflow

Author: Jared G. Galloway
*/

// Using DSL-2
nextflow.enable.dsl=2

// Import a subworkflow to run the BEER enrichment analysis
// https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeR.html
include { edgeR_enrichment } from './edgeR.nf'

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
    #!/usr/bin/env python3

    from phippery.normalize import counts_per_million
    from phippery.utils import *
    
    ds = load("$phip_data")
    counts_per_million(ds)
    dump(ds, "cpm.phip")
    """
}

process size_factors {
    input: path phip_data
    output: path "sf.phip"
    shell:
    """
    #!/usr/bin/env python3

    from phippery.normalize import size_factors
    from phippery.utils import *
    
    ds = load("$phip_data")
    size_factors(ds)
    dump(ds, "sf.phip") 
    """
}

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
    #!/usr/bin/env python3

    from phippery.normalize import enrichment
    from phippery.utils import *
    
    ds = load("$phip_data")
    lib_ds = ds_query(ds, "control_status == 'library'")
    enrichment(ds, lib_ds, data_table="cpm")
    dump(ds, "fold_enr.phip") 
    """
}


process fit_predict_zscore {
    input: path phip_data
    output: path "fit_predict_zscore.phip"
    when: params.run_zscore_fit_predict || params.summarize_by_organism
    shell:
    """
    fit-predict-zscore.py \
        -ds ${phip_data} \
        -o fit_predict_zscore.phip 
    """
}


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

    if( params.run_edgeR_save_rds )
        dataset | edgeR_enrichment | set { edgeR_ch }

    // run some optional statistics which
    // depend on certain annotations
    cpm_fold_enrichment(counts_per_million.out) | set { cpm_fold_enr_ch }
    fit_predict_zscore(counts_per_million.out) | set { fit_pred_zscore_ch }

    // collect all the datasets statistics and merge
    auto_stats_ch.concat(
        cpm_fold_enr_ch,
        fit_pred_zscore_ch,
        edgeR_ch
    ) | collect | merge_binary_datasets

    emit:
    merge_binary_datasets.out
} 
