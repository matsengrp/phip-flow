"""
@Author: Jared Galloway
@date: 5/26/2020

phip-seq data set simulator.

The idea here being the simulation 
gives us an expected outcome from the pipeline,
this is purely for testing purposes at the moment.
"""

import numpy as np
import pandas as pd
import os
from pathlib import Path

nt = ["A","T","C","G"]  


def create_peptide_metadata(
    n_peptides=10,
    tile_length = 93,
    adapt_3_length = 18,
    adapt_5_length = 16,
):
    """
    This function will take in some number of peptides,
    and a tile length to generate random peptides 
    we would expect in a library.
    """

    adapter_3 = ''.join(np.random.choice(nt, adapt_3_length)).lower()
    adapter_5 = ''.join(np.random.choice(nt, adapt_5_length)).lower()
    oligos = []
    for ID in range(n_peptides):
        oligo = ''.join(np.random.choice(nt, tile_length))
        oligo_w_adapt = adapter_3 + oligo + adapter_5
        oligos.append(oligo_w_adapt)

    return oligos, adapter_3, adapter_5


def generate_reads(
    counts,
    samples,
    reference,
    config,
    n_mismatches = 0, 
    replicate_number = 1,
    read_length = 125,
    peptide_on_5_prime = True
):

    sample_fqfp = {}
    for idx, sample in samples.iterrows():
        sample_fqfp[idx] = open(
            os.path.join(
                config[f"experiments"][f"{sample.experiment}"],
                f"{sample.fastq_pattern}".replace("*", f"{replicate_number}")
                ),
            "a"
            )
 
    for pID, peptide in reference.iterrows():
        for sID, sample in samples.iterrows():
            for read_idx in range(counts[pID][sID]):

                # right now we assume the adapters are lower case.
                oligo = generate_mismatch(
                    sequence = "".join(
                        [char for char in peptide["Oligo"] if char.isupper()]
                    ),
                    num_mm = n_mismatches
                )
                
                n_nt_filler = read_length - len(oligo)
                if n_nt_filler > 0:
                    filler = ''.join(np.random.choice(nt, n_nt_filler))
                    split = 0 if peptide_on_5_prime else np.random.choice(range(n_nt_filler))
                    oligo = filler[:split] + oligo + filler[split:]
                qscore = "F" * read_length 
                sample_fqfp[sID].write(f"@\n{oligo}\n+\n{qscore}\n")


def generate_mismatch(sequence, num_mm):
    """
    Take a sequence of nucleotides and mutate (strictly to diff base)
    This should be used to emulate n mismatches during alignment.
    """

    nts = set(["A","T","C","G"])
    mismatch_sequence = ""
    mismatch_indices = np.random.choice(range(len(sequence)), num_mm, replace=False)
    for idx, nt in enumerate(sequence):
        if idx not in mismatch_indices:
            mismatch_sequence += nt
        else:
            mismatch_sequence += np.random.choice(list(nts - set(nt)))
    return mismatch_sequence


def generate_directory_structure(config_dict):
    """
    Take in a dictionary describing the directory structure
    for a simulated dataset and create it
    """

    required = [
        "output_dir",
        "references",
        "experiments",
        "samples"
    ]
    
    assert np.all([req in config_dict for req in required])
    
    base_dir = config_dict["output_dir"]
    for req in required:
        Path(os.path.join(base_dir,req)).mkdir(parents=True, exist_ok=True)
    for experiment in config_dict["experiments"]:
        Path(
            os.path.dirname(config_dict["experiments"][experiment])
        ).mkdir(parents=True, exist_ok=True)
    for experiment in config_dict["references"]:
        Path(
            os.path.dirname(config_dict["references"][experiment])
        ).mkdir(parents=True, exist_ok=True)





