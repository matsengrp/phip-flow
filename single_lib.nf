/*
main.nf

A basic pipeline for creating counts matrices using bowtie alignment 
and the phippery tool

Jared Galloway: 5/26/2020
*/

import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()
params.configFileJS = params.params_file
String configJSON = new File("${params.configFileJS}").text
def config = jsonSlurper.parseText(configJSON)
assert config instanceof Map

peptide_reference_ch = Channel.fromPath(config["peptide_metadata"])

// CONVERT PEPTsample_idE METADATA TO FASTA
process generate_fasta_reference {

    label 'single_thread_large_mem'

    input: file("pep_ref") from peptide_reference_ch

    output:
        set( 
            val(ref_name), 
            file("${ref_name}.fasta") 
        ) into pep_channel_fasta
    
    exec:
        ref_name = config["reference_name"]

    shell:    
        """
        phippery peptide-md-to-fasta -d ${pep_ref} -o ${ref_name}.fasta
        """
}

//peptide_fasta_ch.subscribe{println it}

// GENERATE INDEX
process generate_index {
 
    label 'multithread'
    publishDir "${params.data_dir}/reference"

    //input: file("pep_fasta") from peptide_fasta_ch

    //output: file("peptide_index") into peptide_index_ch
    input:
        set( 
            val(ref_name), 
            file("pep_fasta") 
        ) from pep_channel_fasta

    output:
        set(
            val(ref_name), 
            file("${ref_name}_index") 
        ) into pep_channel_index

    shell:    
        """
        mkdir ${ref_name}_index
        bowtie2-build --threads 4 \
        ${pep_fasta} ${ref_name}_index/${ref_name}
        """
}

//peptide_index_ch.subscribe{println it}

// CREATE SAMPLE CHANNEL
Channel
    .fromPath(config["sample_metadata"])
    .splitCsv(header:true)
    .map{ row -> 
        tuple(
            row.reference,
            row.sample_id,
            new File("${params.data_dir}/NGS/${row.seq_dir}/${row.fastq_filename}"),
        ) 
    }
    .set { samples_ch }
    //.subscribe{println it}

//peptide_index_ch.first().subscribe{println it}


index_sample_ch = pep_channel_index
    .cross(samples_ch)
    .map{ ref, sample ->
        tuple(
            //sample[0],
            //sample[1],
            //fileref[0]
            sample[1],
            ref[0],
            file(ref[1]), // TODO, this is already a file.
            file(sample[2]),
        )
    }

//index_sample_ch.subscribe{println it}

//samples_ch.subscribe{println it}
//peptide_index_v_ch = peptide_index_ch.first()
//index_value_channel = Channel.value(peptide_index_ch.subscribe)

// ALIGN ALL SAMPLES TO THE REFERENCE
process short_read_alignment {

    label 'multithread'
    publishDir "${params.data_dir}/alignments"

    input:
        set( 
            val(sample_id), 
            val(ref_name),
            file(index),
            file(respective_replicate_path),
        ) from index_sample_ch

    output:
        set(
            val(sample_id),
            file("${sample_id}.sam")
        ) into aligned_reads_sam

    exec:
    read_length = config["read_length"]
    tile_length = config["tile_length"]
    trim = read_length - tile_length
    trim = Math.max(0, trim)
    num_mm = config["num_mm"]
    stream_func = config["fastq_stream_func"]

    shell:
        """
        ${stream_func} ${respective_replicate_path} |
        bowtie2 -N ${num_mm} -x ${index}/${ref_name} --trim3 ${trim} --local - > ${sample_id}.sam
        """
            //"""
            //${stream_func} ${respective_replicate_path} |
            //bowtie2 -N ${num_mm} -x ${peptide_index}/index --trim3 ${trim} --local - > ${sample_id}.sam
            //"""
}


//
//aligned_reads_sam.subscribe{println it}
//
//
aligned_reads_sam.into{aligned_reads_for_counts; aligned_reads_for_stats}

//BOWTIE1
process sam_to_stats {

    publishDir "$params.data_dir/stats/"
    label 'multithread'

    input:
        set(
            val(sample_id),
            file(sam_file)
        ) from aligned_reads_for_stats
    
    output:
        file("${sample_id}.txt") into alignment_stats_ch
    
    shell:
        """
        samtools stats ${sam_file} | grep ^SN | cut -f 2- | 
        sed '1p;7p;22p;25p;d' > ${sample_id}.txt
        """ 
}

//alignment_stats_ch.subscribe{println it}

process sam_to_counts {
    
    label 'multithread'
    publishDir "$params.data_dir/counts/"

    input:
        set(
            val(sample_id),
            file(sam_file)
        ) from aligned_reads_for_counts

    output:
        file("${sample_id}.tsv") into counts_ch

    script:
        """
        samtools view -u -@ 4 ${sam_file} | \
        samtools sort -@ 4 - > ${sample_id}.bam
        samtools sort -@ 4 ${sample_id}.bam -o ${sample_id}.sorted 
        mv ${sample_id}.sorted ${sample_id}.bam
        samtools index -b ${sample_id}.bam
        samtools idxstats ${sample_id}.bam | \
        cut -f 1,3 | sed "/^*/d" > ${sample_id}.tsv
        """
}

grouped_counts = counts_ch
    .collect()
    .map{ tsv ->
        tuple(
            file(config["peptide_metadata"]),
            file(config["sample_metadata"]),
            tsv
        )
    }

grouped_stats = alignment_stats_ch.collect()

//TODO not sure how to fix this but "-resume" skips this even 
// if the params inside the configuration have changed.    
process collect_phip_data {
    
    publishDir "$params.data_dir/phip_data/", mode: 'copy'
    label 'single_thread_large_mem'

    input:
        set (
            file(pep_meta),
            file(sam_meta),
            file(all_counts_files)    
        ) from grouped_counts 
        file(all_alignment_stats) from grouped_stats

    output:
        file "${prefix}.phip" into phip_data_ch

    exec:
        prefix = config["dataset_prefix"]

    script:
        """
        phippery collect-phip-data -s_meta ${sam_meta} -p_meta ${pep_meta} \
        -o ${prefix}.phip ${all_counts_files}
        phippery add-stats -ds ${prefix}.phip \
        -o ${prefix}.phip ${all_alignment_stats}
        """ 
}

phip_data_ch.subscribe{println it}


