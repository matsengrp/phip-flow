"""
Generate a single dataset, dubbed test_ones, 
 - no technical replicates, 
 - no duplicate peptide oligos
 - every sample with a read which have a full
    match with each peptide in the library,
    aka the counts matrix should be all ones once
    collected through the pipeline.
"""

import numpy as np
import pandas as pd
from phip_simulator_utils import *
import sys
import json

if __name__ == "__main__":
    
    np.random.seed(23)

    # Custom generate some random samples
    config = json.load(open("simulate_tiny_config.json","r"))

    generate_directory_structure(config)

    fastq_pattern = [f"tiny-sample-*-{i}" for i in range(2)] 
    reference = [f"refa" for _ in range(2)]
    experiment = [f"expa" for _ in range(2)]

    columns = ["ID","fastq_pattern","reference","experiment"]
    df = pd.DataFrame(
        zip(range(len(reference)), fastq_pattern, reference, experiment), columns=columns
    )
    sample_metadata = df.set_index("ID")
    sample_metadata.to_csv(config["samples"], index=True)

    refa, a3, a5 = create_peptide_metadata(
        n_peptides = 2,
        tile_length = 5,
        adapt_3_length = 2,
        adapt_5_length = 2
    )
     
    df = pd.DataFrame(
        zip(range(len(refa)), refa), columns = ["ID", "Oligo"]
    )
    peptide_metadata_refa = df.set_index("ID")
    peptide_metadata_refa.to_csv(config["references"]["refa"], index=True)

    counts = np.ones([len(sample_metadata),len(peptide_metadata_refa)]).astype(int)

    generate_reads(
        counts, 
        sample_metadata, 
        peptide_metadata_refa, 
        config,
        replicate_number = 1,
        read_length = 7)

    

 








