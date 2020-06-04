# phippery-experiments

This repo contains the code for the 
[PhIP-Seq](https://www.nature.com/articles/s41596-018-0025-6), 
analysis pipeline.
Implimented in 
[Nextflow](https://www.nextflow.io/docs/latest/channel.html).

this pipeline requires
(1) a csv sample metadata file specifying demultiplexed NGS fastq files for each sample,
(2) a (number of) csv peptide metadata file(s) specifying the oligo sequence 
for each peptide in the library prepared for Immuno-Precitation (IP).
(3) Finally, the user provides a configuarion file in the form a
[JSON]() 
file to specify file paths to metadata and sequencing files
relative to the working directory of pipeline execution.
We build the index reference and use
[Bowtie]() 
for short read alignment samples to their respective peptide reference library.
The pipeline then merges all alignemnt hit counts into a coherent counts matrix, _M_. 
Concretely, If sample $j$ contains N aligned hits with peptide $i$, then

_M_[i][j] = _N_

This matrix and associated metadata on both axis
can then be used for various statistical queries
or dumped to csv for third-party analysis. 

It's important to note that each sample is associated with a peptide library
for which the IP experiment was run. The respective library for each sample 
in the sample metadata file is defined by a single peptide 
should


## Sample Metadata

A simple csv file containing some metadata for each sample.
Each sample (row) is defined by:

 1. `ID` <int> - A unique identifier to forever tie this sample to it's
    counts for each peptide and metadata for downstream analysis.

 2. `fastq_pattern` <str> - the regex string pattern for the basename of
    of the sample fastq filename.

 3. `experiment` <str> - the "name" of the experiment this sample was run with,
    the config file will tie the this name with a relative path to look for 
    the respective sample filename.

 4. `reference` <str> - the "name" of the reference library this sample should
    be aligned against. 

these are the _required_ fields, but other metadata you would like for downstream
analysis should be tied in here, too. and example might look like:

```
ID,fastq_pattern,experiment,reference,Sample_type,Notes
0,sample-*-0.fastq,lib_ones,refa,_,_
1,sample-*-1.fastq,lib_ones,refa,_,_
2,sample-*-2.fastq,lib_ones,refa,_,_
3,sample-*-3.fastq,lib_ones,refa,_,_
4,sample-*-4.fastq,lib_ones,refa,_,_
5,sample-*-5.fastq,lib_ones,refa,_,_
6,sample-*-6.fastq,lib_ones,refa,_,_
7,sample-*-7.fastq,lib_ones,refa,_,_
8,sample-*-8.fastq,lib_ones,refa,_,_
```

## Peptide Metadata

Another simple csv containing some metadata for each peptide.
Each peptide (row) is defined by:

 1. `ID` <int> - A unique identifier to forever tie this peptide to it's
    counts for each peptide and metadata for downstream analysis.

 2. `Oligo` <str> - the oligo nucleotide sequence defining the expressed
    protein. All Oligos in a peptide_metadata consititute a single library
    used to create the Index. Currently, we only support uppercase representing
    the oligo, every lowercase charictar will be considered an adapter sequence,
    and left out of the Index each sample is aligned against.

these are the _required_ fields, but other metadata you would like for downstream
analysis should be tied in here, too. and example might look like:

```
ID,Oligo
0,gcatcagtaggctgcgtaGGGATTAGGCGGACCTCCATGAATACCGCCATCACAACGCGACCCTGGCTAGCGGCGTTCACGATCAAAGTTACTTTAGTCATGGCTCCATACtcgttaatatgcctgt
1,gcatcagtaggctgcgtaTGTAGGCAAGGAGCAACACTTCTTCTTTGAACTAAGGCTCGCAGAAGTCCCCCATTCTAGCAGGCCGTGCGATCGGGACCGTCGCTTTATTTCtcgttaatatgcctgt
2,gcatcagtaggctgcgtaGAGAATGGGCCAGGAATGATCTACTGTCCTCAATCTTAATAGCATTTGCACTCACTAGGTAAATTCTAAAAATAACTTAATGCGAATTATGCGtcgttaatatgcctgt
3,gcatcagtaggctgcgtaCGTGTCAAAAACTGCGTATTTACGAAGAGATGGTAGAATGGCGGATGTTAAGATAAGACACGGGGCAGGTTGAATTCCATAAAGTTAGTGGAAtcgttaatatgcctgt
4,gcatcagtaggctgcgtaTTTCAGATCCTACCATTTGTGTCCTTAAACGGTCAGAACGTACGAGAGTAGTATGGGGGTTAAGTGTAAGCAAGATCTGACTTGGCGCATGTCtcgttaatatgcctgt
5,gcatcagtaggctgcgtaCCGAGTTCGTATTTTTACAAATCCCGGACTGCATCGTCCTTTCATGTAGCACGGGCCCTGTGTCAGACGCACGATTTCTCCTAGAATTGCTCTtcgttaatatgcctgt
6,gcatcagtaggctgcgtaTATTTAATGAGTGTGAGGCAAAGTTGTTCGGCTCTAGCAAAAGGACGACAAATGAACTAGCCGGAGAACAGCAGTAGTTAAAAGTTATAAGAAtcgttaatatgcctgt
7,gcatcagtaggctgcgtaTTTACGCTCAGCAAGCGTAGCTAGCATGCGTCTTAATGATTCACAACTTTCCTTTATGCATGAACATTCTCTGTCGCTTGGGGGGATGTACTCtcgttaatatgcctgt
8,gcatcagtaggctgcgtaTCAAACAGGTTACGACACAAAGAACGCCAAGTATCTCCGAATCGTACAATCGTGTAGATTTGTTGAGATAGAGTTAACGTAGAGCGCAATTCAtcgttaatatgcctgt
```

## Configuration file.

To specify the relative filepaths for metadata and sequencing files
the user creates and specifies a JSON configuration file. 
The configuration file will specify _one_ sample metadata filepath
and any number of peptide metadata reference libraries associated by
library name. Finally, the file must map an _experiment_ name to a relative
file path we can expect to find the sample fastq files.

```JSON
{
    "output_dir" : "simulations/simulate_ones/",

    "samples" : "simulations/simulate_ones/samples/sample_metadata.csv",

    "experiments" : {
        "expa" : "simulations/simulate_ones/experiments/expa/some/path/",
        "expb" : "simulations/simulate_ones/experiments/expb/some/path/"
    },

    "references" : {
        "refa" : "simulations/simulate_ones/references/peptide_metadata_a.csv",
        "refb" : "simulations/simulate_ones/references/peptide_metadata_b.csv"
    }
}
```

Notice that that _experiments_ and _references_ are both collections themselves. 
This is because the sample metadata is expected to have columns specifying
the associated _experiment_ ("expa", "expb") and _reference_ ("refa", "refb") 
for each sample provided. 





