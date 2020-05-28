# phippery-experiments

This repo contains the code for the PhIP-Seq, analysis pipeline.
Implimented in 
[Nextflow](TODO),
this pipeline takes sample library metadata and sample metadata
in the form of csv files (TODO tsv, too) where the first column
of each file is a unique identifier for the specific sample or peptide 
of interest. The pipeline then reads in associated, 
demiplexed fastq sequencing files for each sample
and aligns them to their respective peptide reference library before 
merging all counts into a coherent counts matrix, $M$. Concretely,
If sample $j$ contains N aligned hits with peptide $i$, then

```
$$M_{i}{j} = N$$
```

This matrix and associated metadata on both axis
can then be used for various statistical queries
or dumped to csv for third-party analysis. 

The last thing required by the pipeline is a JSON
formatted config file specifying the location of required data.

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
Each peptide (row) is defined by

TODO




