# phip-flow

This repo contains the code for our 
[PhIP-Seq](https://www.nature.com/articles/s41596-018-0025-6) 
analysis pipeline.
Implemented in 
[Nextflow](https://www.nextflow.io/docs/latest/channel.html).

this pipeline requires

1. a csv sample metadata file specifying demultiplexed NGS fastq files for each sample,
2. a (possibly number of) csv peptide metadata file(s) specifying the oligo sequence 
for each peptide in the library prepared for Immuno-Precipitation (IP).
3. Finally, the user provides a configuration file in the form a
[JSON]() 
file to specify file paths to metadata and sequencing files
relative to the working directory of pipeline execution.
We build the index reference and use
[Bowtie]() 
for short read alignment samples to their respective peptide reference library.
The pipeline then merges all alignment hit counts into a coherent counts matrix, _M_. 
Concretely, If sample j contains N aligned hits with peptide i, then _M_[i][j] = _N_.

This matrix and associated metadata on both axes
can then be used for various statistical queries
or dumped to csv for third-party analysis. 

Each sample is associated with a peptide library
for which the IP experiment was run. The respective library for each sample
in the sample metadata file is defined by a single reference.
You can mix as many samples and references as you'd like in a single sample
metadata table, but note that the pipeline will create an enrichment
dataset for each reference, each containing only the relevent samples.

## Sample Metadata

A simple csv file containing some metadata for each sample.
Each sample (row) is defined by:

 1. `ID` <int> - A unique identifier to forever tie this sample to it's
    counts for each peptide and metadata for downstream analysis.

 2. `fastq_pattern` <str> - the regex string pattern for the basename of
    of the sample fastq filename.

 3. `seq_dir` <str> - the "name" of the sequencing sir this sample was run with,
    the config file will tie the this name with a relative path to look for 
    the respective sample filename.

 4. `reference` <str> - the "name" of the reference library this sample should
    be aligned against. 

these are the _required_ fields, but other metadata you would like for downstream
analysis should be tied in here, too. and example might look like:

```
ID,fastq_pattern,reference,seq_dir
0,sample-*-0,refa,expa
1,sample-*-1,refa,expa
2,xeno-AE-122-*-R1.3.3.0,refa,expa
3,xeno-AE-122-*-R1.3.3.1,refa,expa
4,xeno-AE-122-*-R1.3.3.2,refa,expa
5,xeno-AE-122-*-R1.3.3.3,refb,expb
6,xeno-AE-122-*-R1.3.3.4,refb,expb
7,johnny-boy-*-0,refb,expb
8,johnny-boy-*-1,refb,expb
9,johnny-boy-*-2,refb,expb
```

## Peptide Metadata

Another simple csv containing some metadata for each peptide.
Each peptide (row) is defined by:

 1. `ID` <int> - A unique identifier to forever tie this peptide to its
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

## Configuration file

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

    "seq_dir" : {
        "expa" : "simulations/simulate_ones/experiments/expa/some/path/",
        "expb" : "simulations/simulate_ones/experiments/expb/some/path/"
    },

    "references" : {
        "refa" : "simulations/simulate_ones/references/peptide_metadata_a.csv",
        "refb" : "simulations/simulate_ones/references/peptide_metadata_b.csv"
    },

    "read_length" : {
        "seq_dir_1" : 125,
    },

    "tile_length" : {
        "cov" : 117
    },

    "num_mm" : 2,

    "counts_matrix_prefix" : "Example-",
    
    "tech_rep_agg_func" : "mean",

    "fastq_stream_func" : "zcat"
}
```

Notice that that _experiments_ and _references_ are both collections themselves. 
This is because the sample metadata is expected to have columns specifying
the associated _experiment_ ("expa", "expb") and _reference_ ("refa", "refb") 
for each sample provided. 

# Example running on gizmo

To run on gizmo, I like to start by creating a directory for the inputs, config files
and where the outputs will be stored. For example such a directory might look like this:

```
(base) ➜  pipeline-run-07-17-20 (master) tree
.
├── config.json             # the config file for nextflow to find inputs
├── nextflow.gizmo.config   # the cluster config file specifying partitions, containers etc.
├── peptide_metadata.csv    # peptide metadata for the run
├── run_gizmo.sh            # a shell script with the nextflow run command
└── sample_metadata.csv     # sample metadata for the run

0 directories, 5 files 
```

With this directory, the `config.json` file might look something like:

```json
{
    "output_dir" : "./",

    "samples" : "sample_metadata.csv",

    "seq_dir" : {
        "seq_dir_1" : "/path/to/seq_dir/",
    },

    "references" : {
        "cov" : "peptide_metadata.csv"
    },

    "read_length" : {
        "seq_dir_1" : 125,
    },

    "tile_length" : {
        "cov" : 117
    },

    "num_mm" : 2,

    "counts_matrix_prefix" : "07-17-20-",
    
    "tech_rep_agg_func" : "mean",

    "fastq_stream_func" : "zcat"
}
```

The shell script for running nextflow might look something like this:
```
#!/bin/bash

set -e
source /app/lmod/lmod/init/profile

module load nextflow
module load Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

/usr/bin/time nextflow \
    run \
    ../../phip-flow/phip-flow/main.nf \     # path to where you cloned this repo containing the main.nf script.
    -c nextflow.gizmo.config \
    --params_file config.json \
    -with-report output/nextflow_report.html \
    -work-dir /fh/scratch/delete30/some/temp/path/to/store/outputs/ \ # the path where all intermediate files are stored
    -ansi-log false \
```


Next, on rhino, request an interactive node where you will run the pipeline.

```
(base) ➜  pipeline-run-07-17-20 (master) kinit
Password for jgallowa@FHCRC.ORG: <insert password>
(base) ➜  pipeline-run-07-17-20 (master) grabnode
How many CPUs/cores would you like to grab on the node? [1-36] 12
How much memory (GB) would you like to grab? [240] 16
Please enter the max number of days you would like to grab this node: [1-7] 1
Do you need a GPU ? [y/N]N
```

Once allocated, navigate to the directory where you put the above files and simply

```
./run_gizmo.sh
```

I reccomend doing this in a `tmux` shell incase your connection to the server breaks before the job finishes.



