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
    config = json.load(open("simulate_ones_single_ref.json","r"))

    generate_directory_structure(config)
    
    fastq_pattern = [f"xeno-AE-122-*-R1.3.3.{i}" for i in range(5)]
    reference = [f"refa" for _ in range(5)] 
    experiment = [f"expa" for _ in range(5)]

    columns = ["ID","fastq_pattern","reference","experiment"]
    df = pd.DataFrame(
        zip(range(len(reference)), fastq_pattern, reference, experiment), columns=columns
    )
    sample_metadata = df.set_index("ID")
    sample_metadata.to_csv(config["samples"], index=True)

    refa, a3, a5 = create_peptide_metadata(tile_length = 117)
    df = pd.DataFrame(
        zip(range(len(refa)), refa), columns = ["ID", "Oligo"]
    )
    peptide_metadata_refa = df.set_index("ID")
    peptide_metadata_refa.to_csv(config["references"]["refa"], index=True)

    counts = np.ones([len(peptide_metadata_refa),len(sample_metadata)]).astype(int)
    print(counts)

    generate_reads(
        counts,
        sample_metadata, 
        peptide_metadata_refa, 
        config,
        read_length=125,
        replicate_number=1)

    generate_reads(
        counts,
        sample_metadata, 
        peptide_metadata_refa, 
        config,
        read_length=125,
        replicate_number=2)
