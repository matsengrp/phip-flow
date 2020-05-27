log.info('main.nf is running')

import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()

params.configFileJS = "config.json"
String configJSON = new File("${params.configFileJS}").text
def config = jsonSlurper.parseText(configJSON)
//log.info("${config}")
//log.info("${config.references.refa}")
assert config instanceof Map

//config.references.each { entry ->
// println entry.key
// println entry.value
//}

/*
get_values = { ref_key, ref_name_v ->
    ref_name_v
}
x = config["references"].values()
println x
println params.references

*/
//peptideMetadataChannel = Channel.fromPath(params.references.values())

references = config["references"].values()
println references.getClass()
pep_md_fp = new ArrayList<String>(config["references"].values());
pep_channel = Channel.fromPath(pep_md_fp)
//pep_channel.subscribe { println "value: $it" }
//ref = Channel.from("${config.references.refa}")

// you should pss ref name and file.
process create_fasta_reference {

    
    container 'phippery:latest'    
    
    input:
    file pep_ref from pep_channel

    output:
    file "${pep_ref}.fasta" into pep_channel_fasta
    
    shell:    
    """
    phippery peptide-md-to-fasta -d ${pep_ref} -o ${pep_ref}.fasta
    """

}

// this is going to take in all fastq files from samples.
// input channel for fastq should be all samples split by
// reference library.
process create_index {

    
    container 'biocontainers/bowtie:v1.2.2dfsg-4-deb_cv1'    

    input:
    file pep_fasta from pep_channel_fasta

    //output:
    //pep_channel_Index
    
    shell:    
    """
    bowtie-build
    """

}























pep_channel_fasta.println()
