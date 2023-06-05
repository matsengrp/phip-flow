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
    publishDir "$params.results/rds_data/", mode: 'copy', overwrite: true
    input:
    tuple path(phip_data), path(phip_data_csvs)
    output:
    tuple path(phip_data), path("edgeR*.csv"), path("*.rds")
    shell:    
    """
    run_edgeR.Rscript
    mv PhIPData.rds ${params.dataset_prefix}.rds
    """
}

// APPEND EDGER RESULTS INTO XARRAY DATASET
process append_edgeR_to_xarray {
    input:
    tuple path(phip_data), path(edgeR_csvs), path(rds_data)
    output:
    path "edgeR.phip" 
    shell:
    """
    #!/usr/bin/env python3

    import glob
    from phippery.utils import *
    import pandas as pd

    ds = load("$phip_data")
    for edgeR_csv in glob.glob("*.csv"):
        df = pd.read_csv(edgeR_csv, index_col=0)
        table_name = edgeR_csv.split(".")[0]
        add_enrichment_layer_from_array(
            ds, df.values, new_table_name=table_name
        )

    dump(ds, "edgeR.phip") 
    """
}

workflow edgeR_enrichment {
    take:
        ds
    main:

    ds | to_csv | run_edgeR | append_edgeR_to_xarray

    emit:
        append_edgeR_to_xarray.out

}


