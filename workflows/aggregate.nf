#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


// Get the list of all samples so that aggregate_organisms can be sharded
process split_samples {

    input:
        // All output data in wide format (CSV)
        path "*"

    output: path "sample_list"
    when: params.summarize_by_organism
    shell:
    template "split_samples.py"
}

process aggregate_organisms {
    tag "${sample_id}"
    cpus 1
    memory "4.GB"
    input:
        // All output data in wide format (CSV)
        tuple path("*"), val(sample_id)
        // Any public epitopes defined in CSV format
        path public_epitopes_csv
    output: path "*.csv.gz"
    when: params.summarize_by_organism
    shell:
    template "aggregate_organisms.py"
}

process join_organisms {
    publishDir "$params.results/aggregated_data/", mode: 'copy', overwrite: true
    input: path "input/"
    output: path "*.csv.gz"
    when: params.summarize_by_organism
    shell:
    template 'join_organisms.py'
}

workflow AGG {
    take:
        dump_binary
        dump_wide_csv
        dump_tall_csv
    main:

    // Get the list of all samples
    split_samples(dump_wide_csv)

    aggregate_organisms(
        dump_wide_csv
            .toSortedList()
            .combine(
                split_samples
                    .out
                    .splitText(){it.replace("\n", "")}
            ),
        file("${params.public_epitopes_csv}")
    )

    join_organisms(
        aggregate_organisms
            .out
            .flatten()
            .toSortedList()
    )
}


