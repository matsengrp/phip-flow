/*
main.nf

A basic pipeline for creating counts matrices using bowtie alignment 
and the phippery tool

Jared Galloway: 5/26/2020
*/


import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()
params.configFileJS = "config.json"
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

index_sample_ch = pep_channel_index
    .cross(samples_ch)
    .map{ ref, sample ->
        tuple(
            sample[1],
            ref[0],
            file(ref[1]),
            file(sample[2])
        )
    }
    //.subscribe {
    //    println it
    //    it.each { arg ->
    //        println arg.getClass()
    //    }
    //}

process short_read_alignment {

    publishDir "$config.output_dir/references/"
    container 'quay.io/biocontainers/bowtie:1.2.2--py36h2d50403_1'    
    echo true

    input:
        set( 
            val(ID), 
            val(ref_name), 
            file(index), 
            file(technical_replicates)
        ) from index_sample_ch 

    shell:
    """
    bowtie ${index}/${ref_name} ${technical_replicates.toString()} --sam
    """    
}












