/*
main.nf

A basic pipeline for creating counts matrices using bowtie alignment 
and the phippery tool

Jared Galloway: 5/26/2020
*/


log.info('main.nf is running')

// parse the config file giving us 
// details about the 
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()
params.configFileJS = "config.json"
String configJSON = new File("${params.configFileJS}").text
def config = jsonSlurper.parseText(configJSON)
assert config instanceof Map

references = config["references"].values()
println references.getClass()
pep_md_fp = new ArrayList<String>(config["references"].values());
pep_channel = Channel.fromPath(pep_md_fp)
//pep_channel.subscribe { println "value: $it" }



// TODO: you should pss ref name and file? for index naming and whatnot
// TODO: should we have an option to publish the directory somewhere for int. files?
process create_fasta_reference {

    // TODO publish this on quay
    container 'phippery:latest'    

    input:
    file pep_ref from pep_channel

    output:
    file "${pep_ref}.fasta" into pep_channel_fasta

    // TODO let's strip the .csv fe off of the fasta file    
    shell:    
    """
    phippery peptide-md-to-fasta -d ${pep_ref} -o ${pep_ref}.fasta
    """
}


// this is going to take in all fastq files from samples.
// input channel for fastq should be all samples split by
// reference library.
process create_index {
 
    container 'quay.io/biocontainers/bowtie:1.2.2--py36h2d50403_1'    

    input:
    file pep_fasta from pep_channel_fasta.flatten()

    output:
    file "Index" into pep_channel_Index
   
    // TODO  
    shell:    
    """
    mkdir Index
    bowtie-build ${pep_fasta} Index/${pep_fasta}
    """

}

pep_channel_Index.subscribe {println it}



