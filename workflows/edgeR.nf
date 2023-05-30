// Run external statistical analysis tools


// EXTRACT WIDE CSV
process to_csv {
    input: path phip_data
    output: path "*.csv"
    shell:
    """
    phippery to-wide-csv -o dataset $phip_data
    """
}

// RUN BEER
process run_edgeR {
    publishDir "$params.results/rds_data/", mode: 'copy', overwrite: true
    input:
    path "*"
    output:
    path "*PhIPData.rds"
    shell:    
    """
    run_edgeR.Rscript
    mv PhIPData.rds ${params.dataset_prefix}-PhIPData.rds
    """
}

workflow edgeR_enrichment {
    take:
        ds
    main:

    ds | to_csv | run_edgeR

    emit:
        run_edgeR.out

}


