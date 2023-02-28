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
process run_beer {
    input:
    path "*"
    output:
    path "beer_results.tsv"
    shell:    
    """
    run_beer.Rscript
    """
}

workflow beer_enrichment {
    take:
        ds
    main:

    ds | to_csv | run_beer

    emit:
        run_beer.out

}


