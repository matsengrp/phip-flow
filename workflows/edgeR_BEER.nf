// Run external statistical analysis tools


// EXTRACT WIDE CSV
process to_csv {
    input: path phip_data
    output: 
    tuple path(phip_data), path("*.csv")
    shell:
    """
    phippery to-wide-csv -o dataset $phip_data
    """
}

// RUN BEER
process run_edgeR {
    // publishDir "$params.results/rds_data/", mode: 'copy', overwrite: true
    input:
    tuple path(phip_data), path(phip_data_csvs)
    output:
    tuple path(phip_data), path("edgeR*.csv"), path("PhIPData.rds"), val("edgeR")
    shell:    
    """
    run_edgeR.Rscript ${params.edgeR_threshold}
    """
}
//mv PhIPData.rds ${params.dataset_prefix}.rds

process run_BEER {
    // publishDir "$params.results/rds_data/", mode: 'copy', overwrite: true
    input:
    tuple path(phip_data), path("*"), path(edgeR_rds), val(method)
    output:
    tuple path(phip_data), path("beer*.csv"), path("PhIPData.rds"), val("BEER")
    shell:    
    """
    run_BEER.Rscript
    """

}

process publish_rds {
    publishDir "$params.results/rds_data/", mode: 'copy', overwrite: true
    input:
    tuple path(phip_data), path(csvs), path(rds_data), val(method)
    output:
    path rds_data
    """
    echo publishing $rds_data 
    """
}

// APPEND EDGER RESULTS INTO XARRAY DATASET
process append_assay_csvs_to_xarray {
    input:
    tuple path(phip_data), path(csvs), path(rds_data), val(method)
    output:
    path "${method}.phip"
    shell:
    """
    #!/usr/bin/env python3

    import glob
    from phippery.utils import *
    import pandas as pd

    ds = load("$phip_data")
    for csv in glob.glob("*.csv"):
        df = pd.read_csv(csv, index_col=0)
        table_name = csv.split(".")[0]
        add_enrichment_layer_from_array(
            ds, df.values, new_table_name=table_name
        )

    dump(ds, "${method}.phip") 
    """
}

workflow edgeR_BEER_workflows {
    take:
        ds
    main:

    if ( params.run_BEER )
        ds | to_csv \
            | run_edgeR \
            | run_BEER \
            | (append_assay_csvs_to_xarray & publish_rds)
    else
        ds | to_csv \
            | run_edgeR \
            | (append_assay_csvs_to_xarray & publish_rds)

    emit:
        append_assay_csvs_to_xarray.out

}


