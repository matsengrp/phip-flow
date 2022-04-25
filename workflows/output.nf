#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process dump_tall_csv {
    publishDir "$params.results/tall_data/", mode: 'copy', overwrite: true
    input: file phip_data
    output: file "*.csv"
    when: params.output_tall_csv
    shell:
    """
    phippery to-tall-csv -o ${params.dataset_prefix}-tall.csv $phip_data 
    """
}

process dump_wide_csv {
    publishDir "$params.results/wide_data/", mode: 'copy', overwrite: true
    input: path phip_data
    output: path "*.csv"
    when: params.output_wide_csv
    shell:
    """
    phippery to-wide-csv -o $params.dataset_prefix $phip_data
    """
}

process dump_binary {
    publishDir "$params.results/pickle_data/", mode: 'copy', overwrite: true
    input: file phip_data
    output: file "${params.dataset_prefix}.phip"
    when: params.output_pickle_xarray
    shell:
    """
    cp ${phip_data} ${params.dataset_prefix}.phip
    """
}

process aggregate_organisms {
    publishDir "$params.results/aggregated_data/", mode: 'copy', overwrite: true
    input:
        // All output data in wide format (CSV)
        path "*"
        // Any public epitopes defined in CSV format
        path public_epitopes_csv
    output: path "*.csv"
    when: params.summarize_by_organism
    shell:
    template: "aggregate_organisms.py"
}

workflow DSOUT {
    take: dataset
    main:
    dump_binary(dataset)
    dump_wide_csv(dataset)
    dump_tall_csv(dataset)
    aggregate_organisms(
        dump_wide_csv.out,
        file("${params.public_epitopes_csv}")
    )
}


