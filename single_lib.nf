/*
main.nf

A basic pipeline for creating counts matrices using bowtie alignment 
and the phippery tool

Jared Galloway: 5/26/2020
*/


println params.peptide_metadata
peptide_reference_ch = Channel.fromPath("${params.output}/$params.peptide_metadata")


// CONVERT PEPTIDE METADATA TO FASTA
process generate_fasta_reference {

    label 'single_thread_large_mem'
    publishDir "${params.output}/intermediate_file_links/${params.phip_data_dir}/"

    input: file("pep_ref") from peptide_reference_ch

    output: file("peptide.fasta") into pep_channel_fasta

    shell:    
        """
        phippery peptide-md-to-fasta -d ${pep_ref} -o peptide.fasta
        """
}


// GENERATE INDEX
process generate_index {
 
    label 'multithread'
    publishDir "${params.output}/intermediate_file_links/${params.phip_data_dir}/reference"

    input:
        file("pep_fasta") from pep_channel_fasta

    output:
        set(
            val("peptide_ref"), 
            file("peptide_index") 
        ) into pep_channel_index

    shell:    
        """
        mkdir peptide_index
        bowtie2-build --threads 4 \
        ${pep_fasta} peptide_index/peptide
        """
}


// CREATE SAMPLE CHANNEL
Channel
    .fromPath("${params.sample_metadata}")
    .splitCsv(header:true)
    .map{ row -> 
        tuple(
            "peptide_ref",
            row.sample_id,
            new File("${params.output}/NGS/${row.seq_dir}/${row.fastq_filename}"),
        ) 
    }
    .set { samples_ch }

index_sample_ch = pep_channel_index
    .cross(samples_ch)
    .map{ ref, sample ->
        tuple(
            sample[1],          // sample_id
            file(ref[1]),       // index files
            file(sample[2]),    // sample path
        )
    }


// ALIGN ALL SAMPLES TO THE REFERENCE
process short_read_alignment {

    label 'multithread'
    publishDir "${params.output}/intermediate_file_links/${params.phip_data_dir}/alignments"

    input:
        set( 
            val(sample_id), 
            file(index),
            file(respective_replicate_path),
        ) from index_sample_ch

    output:
        set(
            val(sample_id),
            file("${sample_id}.sam")
        ) into aligned_reads_sam

    exec:
        read_length = params.read_length
        tile_length = params.tile_length
        trim = read_length - tile_length
        trim = Math.max(0, trim)
        num_mm = params.num_mm
        stream_func = params.fastq_stream_func

    shell:
        """
        ${stream_func} ${respective_replicate_path} |
        bowtie2 -N ${num_mm} -x ${index}/peptide --trim3 ${trim} - > ${sample_id}.sam
        """
}


// SPLIT CHANNEL FOR COUNTS AND STATS IN PARALLEL
aligned_reads_sam.into{aligned_reads_for_counts; aligned_reads_for_stats}


// ALIGNMENT STATS
process sam_to_stats {

    label 'multithread'
    publishDir "${params.output}/intermediate_file_links/${params.phip_data_dir}/stats"

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


// COMPUTE COUNTS FOR ALL SAMPLES
process sam_to_counts {
    
    publishDir "${params.output}/intermediate_file_links/${params.phip_data_dir}/counts/"
    label 'multithread'

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


// GRAB ALL
grouped_counts = counts_ch
    .collect()
    .map{ tsv ->
        tuple(
            file("${params.peptide_metadata}"),
            file("${params.sample_metadata}"),
            tsv
        )
    }

grouped_stats = alignment_stats_ch.collect()

// COLLECT AND MERGE ALL 
process collect_phip_data {
    
    label 'single_thread_large_mem'
    publishDir "${params.output}/intermediate_file_links/${params.phip_data_dir}/", mode: 'copy'

    input:
        set (
            file(pep_meta),
            file(sam_meta),
            file(all_counts_files)    
        ) from grouped_counts 
        file(all_alignment_stats) from grouped_stats

    output:
        file "${params.dataset_prefix}.phip" into phip_data_ch


    script:
        """
        phippery collect-phip-data -s_meta ${sam_meta} -p_meta ${pep_meta} \
        -o ${params.dataset_prefix}.phip ${all_counts_files}
        phippery add-stats -ds ${params.dataset_prefix}.phip \
        -o ${params.dataset_prefix}.phip ${all_alignment_stats}
        """ 
}

phip_data_ch.subscribe{println it}


