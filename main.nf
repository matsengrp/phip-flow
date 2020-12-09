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

ref_tuple_channel = Channel
    .fromList(
        config["references"].collect { key, val ->
            new Tuple(key, file(val))
        }
    )

process generate_fasta_reference {

    label 'single_thread_large_mem'

    input:
        set( 
            val(ref_name), 
            file("pep_ref") 
        ) from ref_tuple_channel

    output:
        set( 
            val(ref_name), 
            file("${ref_name}.fasta") 
        ) into pep_channel_fasta

    shell:    
    """
    phippery peptide-md-to-fasta -d ${pep_ref} -o ${ref_name}.fasta
    """
}


process generate_index {
 
    label 'multithread'

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
    bowtie-build --threads 4  ${pep_fasta} ${ref_name}_index/${ref_name}
    """
}


Channel
    .fromPath(config["samples"])
    .splitCsv(header:true)
    .map{ row -> 
        tuple(
            row.reference,
            row.sample_id,
            new File("${params.output}/NGS/${row.seq_dir}/${row.fastq_filename}"),
            row.seq_dir
        ) 
    }
    .set { samples_ch }

index_sample_ch = pep_channel_index
    .cross(samples_ch)
    .map{ ref, sample ->
        tuple(
            sample[1],
            ref[0],
            file(ref[1]), // TODO, this is already a file.
            file(sample[2]),
            sample[3] //exp name
        )
    }


process short_read_alignment {

    label 'multithread'

    input:
        set( 
            val(ID), 
            val(ref_name), 
            file(index), 
            file(respective_replicate_path),
            val(experiment_name)
        ) from index_sample_ch 

    output:
        set(
            val(ID),
            val(ref_name),
            file("${ID}.sam")
        ) into aligned_reads_sam

    exec:
    read_length = config["read_length"]["${experiment_name}"] 
    tile_length = config["tile_length"]["${ref_name}"]
    trim = read_length - tile_length
    num_mm = config["num_mm"]
    stream_func = config["fastq_stream_func"]

    shell:
        """
        ${stream_func} ${respective_replicate_path} |
        bowtie --trim3 ${trim} -n ${num_mm} -l ${tile_length} \
        --tryhard --nomaqround --norc --best --sam --quiet \
        ${index}/${ref_name} - \
        > ${ID}.sam
        """
}

aligned_reads_sam.into{aligned_reads_for_counts; aligned_reads_for_stats}

process sam_to_stats {

    publishDir "$params.output/stats/$ref_name"
    label 'multithread'

    input:
        set(
            val(ID),
            val(ref_name),
            file(sam_file)
        ) from aligned_reads_for_stats
    
    output:
        set(
            val(ID),
            val(ref_name),
            file("${ID}.txt")
        ) into alignment_stats 
    
    shell:
        """
        samtools stats ${sam_file} | grep ^SN | cut -f 2- | 
        sed '1p;7p;24p;31p;d' > ${ID}.txt
        """ 
        //cut -f 2 -d\$'\t' | sed '1p;7p;24p;31p;d' > ${ID}.txt

}


process sam_to_counts {
    
    label 'multithread'

    input:
        set(
            val(ID),
            val(ref_name),
            file(sam_file)
        ) from aligned_reads_for_counts

    output:
        set(
            val(ref_name),
            file("${ID}.tsv")
        ) into counts

    script:
        """
        samtools view -u -@ 4 ${sam_file} | \
        samtools sort -@ 4 - > ${ID}.bam
        samtools sort -@ 4 ${ID}.bam -o ${ID}.sorted 
        mv ${ID}.sorted ${ID}.bam
        samtools index -b ${ID}.bam
        samtools idxstats ${ID}.bam | \
        cut -f 1,3 | sed "/^*/d" > ${ID}.tsv
        """
}

grouped_counts = counts
    .groupTuple()
    .map{ ref_name, tsv ->
        tuple(
            ref_name,
            file(config["references"][ref_name]),
            file(config["samples"]),
            tsv
        )
    }

//grouped_counts.subscribe{println it}

//TODO not sure how to fix this but "-resume" skips this even 
// if the params inside the configuration have changed.    
process collect_phip_data {
    
    publishDir "$params.output/phip_data/", mode: 'copy'
    label 'single_thread_large_mem'

    input:
        set (
            val(ref_name),
            file(pep_meta),
            file(sam_meta),
            file(all_counts_files)    
        ) from grouped_counts 

    output:
        file "${prefix}${ref_name}.phip" into phip_data_ch

    exec:
        prefix = config["counts_matrix_prefix"]
        tech_rep_agg_func = config["tech_rep_agg_func"]

    script:
    """
    phippery collect-phip-data -s_meta ${sam_meta} -p_meta ${pep_meta} \
    -o ${prefix}${ref_name}.phip ${all_counts_files}
    """ 
}

phip_data_ch.subscribe{println it}
