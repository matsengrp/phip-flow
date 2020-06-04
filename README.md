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
ID,fastq_pattern,reference,experiment
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

## Testing pipeline with simulation

To test the pipeline and hueristics within the pipeline, it is helpful
to simulate some sequencing files and metadata for which we 
know the expected result (counts table) from the pipeline, given this data.
The file `simulations/phip_simulator_utils.py` provides some useful
functions to simulate a dataset. The user then writes a script 
using these utilities to define reads for
any sample/peptide combination by supplying a matrix of counts for each
n-mismatch profile.

To sun the simulation, you simply need `python3`, `numpy`, and `pandas`.
A conda environment can be created for this using the `environment.yaml`
file in the simulations directory, like so:

```
cd simulations && conda env create -f environment.yaml
conda activate phip-simulation
```

Then, create a custom script and config file using the phip_simulation_utils
provided. Some examples of these can be found in the simulations directory.








