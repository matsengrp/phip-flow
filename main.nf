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

    publishDir "$config.output_dir/references/"
    container 'phippery:latest'    

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
 
    publishDir "$config.output_dir/references/"
    container 'quay.io/biocontainers/bowtie:1.2.2--py36h2d50403_1'    

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
    bowtie-build ${pep_fasta} ${ref_name}_index/${ref_name}
    """
}


Channel
    .fromPath(config["samples"])
    .splitCsv(header:true)
    .map{ row -> 
        tuple(
            row.reference,
            row.ID,
            new File(config["experiments"][row.experiment], row.fastq_pattern)
        ) 
    }
    .set { samples_ch }


// TODO, this could be cleaned up so .flatMap does all the work,
// and we can get rid of .map all together.
index_sample_ch = pep_channel_index
    .cross(samples_ch)
    .map{ ref, sample ->
        tuple(
            sample[1],
            ref[0],
            file(ref[1]), // TODO, this is already a file.
            file(sample[2])
        )
    }
    .flatMap{t ->
        t[3].withIndex().collect{ replicate_path, replicate_idx ->
            tuple(
                t[0], //id
                replicate_idx,
                t[1], //ref name
                t[2], //ref dir
                replicate_path
            )
        }
    }
    //.subscribe { 
        //println "$it"
        //it.each { arg ->
        //    println arg.getClass()
        //}
    //}

process short_read_alignment {

    publishDir "$config.output_dir/alignments/$ref_name/"
    container 'quay.io/biocontainers/bowtie:1.2.2--py36h2d50403_1'
    echo true 

    input:
        set( 
            val(ID), 
            val(replicate_number),
            val(ref_name), 
            file(index), 
            file(respective_replicate_path)
        ) from index_sample_ch 

    output:
        set(
            val(ID),
            val(replicate_number),
            val(ref_name),
            file("${ID}.${replicate_number}.sam")
        ) into aligned_reads_ch

    // TODO obviously we are going to need different alignment
    // schemes for phage-dma, phip-seq, simulation testing.
    // TODO should we be doing local alignment?
    // like, aren't the reads going to be longer than the reference library
    // peptides?
    shell:
    """
    echo blacksheep > ${ID}.${replicate_number}.sam
    echo Unique sample ID: ${ID}
    echo replicate number: ${replicate_number}
    echo reference name: ${ref_name}
    echo reference dir location: ${index}
    """    
    //bowtie -v 2 --tryhard --sam \
    //${index}/${ref_name} ${technical_replicates.toString()} > ${ID}.sam;
}

//process generate_counts {
    
//    publishDir "$config.output_dir/counts/$ref_name"
//    container ''
//}








