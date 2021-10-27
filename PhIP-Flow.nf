#!/usr/bin/env nextflow
/*
PhIP-Flow Pipeline Script

Author: Jared G. Galloway
*/

// Using DSL-2
nextflow.enable.dsl=2

// Initialize parameters used in the workflow
params.dataset_prefix = 'output'
params.phip_data_dir = false
params.peptide_table = false
params.sample_table = false

// CONVERT PEPTIDE METADATA TO FASTA
process generate_fasta_reference {

    label 'single_thread_large_mem'
    //publishDir "${params.phip_data_dir}/"

    input:
        file "pep_ref"

    output:
        file "peptides.fasta"

    shell:
        """
        phipflow peptide-md-to-fasta -d ${pep_ref} -o peptides.fasta
        """
}


// GENERATE INDEX
process generate_index {
 
    //publishDir "${params.phip_data_dir}/reference"
    label 'multithread'

    input:
        file "pep_fasta"

    output:
        tuple val("peptide_ref"), path("peptide_index")

    shell:    
    template "generate_index.sh"

}


// ALIGN ALL SAMPLES TO THE REFERENCE
process short_read_alignment {

    //publishDir "${params.phip_data_dir}/alignments"
    label 'multithread'

    input:
        tuple val(sample_id), path(index), path(respective_replicate_path)

    output:
        tuple val(sample_id), path("${sample_id}.sam")

    shell:
    template "short_read_alignment.sh"

}


// COMPUTE ALIGNMENT STATS FOR ALL STATS
process sam_to_stats {

    //publishDir "${params.phip_data_dir}/stats"
    label 'multithread'

    input:
        tuple val(sample_id), path(sam_file)
    
    output:
        path "${sample_id}.stats"
    
    shell:
    template "sam_to_stats.sh"

}


// COMPUTE COUNTS FOR ALL SAMPLES
process sam_to_counts {
    
    //publishDir "${params.phip_data_dir}/counts/"
    label 'multithread'

    input:
        tuple val(sample_id), file(sam_file)

    output:
        path "${sample_id}.counts"

    script:
    template "sam_to_counts.sh"
}


// COLLECT AND MERGE ALL 
process collect_phip_data {
    
    publishDir "${params.phip_data_dir}/", mode: 'copy', overwrite: true
    label 'single_thread_large_mem'

    input:
        file all_counts_files
        file all_alignment_stats
        file sample_table 
        file peptide_table 

    output:
        path "${params.dataset_prefix}.phip"

    script:
    template "collect_phip_data.sh"
}

// Entrypoint for the main workflow, which is run by default
workflow {

    // Make sure that the user provided a peptide_table parameter
    if (!params.peptide_table){
        log.info"""
        User must provide --peptide_table
        """.stripIndent()
        exit 1
    }

    // Make sure that the user provided a phip_data_dir parameter
    if (!params.phip_data_dir){
        log.info"""
        User must provide --phip_data_dir
        """.stripIndent()
        exit 1
    }

    // Reference the peptide table provided by the user
    peptide_reference_ch = Channel.fromPath("${params.peptide_table}")

    // Convert peptide metadata to FASTA
    generate_fasta_reference(
        peptide_reference_ch
    )

    // Generate an index from that FASTA
    generate_index(
        generate_fasta_reference.out
    )

    // Make sure that the user provided a sample_table parameter
    if (!params.sample_table){
        log.info"""
        User must provide --sample_table
        """.stripIndent()
        exit 1
    }

    // Create a channel with input data provided by the user
    Channel
        .fromPath("${params.sample_table}")
        .splitCsv(header:true)
        .map{ row -> 
            tuple(
                "peptide_ref",
                row.sample_id,
                file("${row.seq_dir}/${row.fastq_filename}")
            ) 
        }
        .set { samples_ch }

    // Align all samples to the reference
    short_read_alignment(
        generate_index.out
            .cross(samples_ch)
            .map{ ref, sample ->
                tuple(
                    sample[1],          // sample_id
                    file(ref[1]),       // index files
                    file(sample[2]),    // sample path
                )
            }
    )

    // Compute counts for all samples
    sam_to_counts(
        short_read_alignment.out
    )

    // Compute alignment stats
    sam_to_stats(
        short_read_alignment.out
    )

    // Collect and merge all
    collect_phip_data(
        sam_to_counts.out.toSortedList(),
        sam_to_stats.out.toSortedList(),
        Channel.fromPath("${params.sample_table}"),
        Channel.fromPath("${params.peptide_table}")
    )

}
