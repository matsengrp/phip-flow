#!/usr/bin/env python3
import os
import pandas as pd

peptide_dtypes = dict(
    peptide=str,
    n_replicates=int,
    EBS=float,
    hit=str,
    sample=str,
    public=bool
)

org_dtypes = dict(
    sample=str,
    organism=str,
    n_hits_all=int,
    n_discordant_all=int,
    max_ebs_all=float,
    mean_ebs_all=float,
    n_hits_public=int,
    n_discordant_public=int,
    max_ebs_public=float,
    mean_ebs_public=float,
    
)

for suffix, dtype_dict in [
    ("peptide.ebs.csv.gz", peptide_dtypes),
    ("organism.summary.csv.gz", org_dtypes)
]:

    pd.concat([
        pd.read_csv(
            os.path.join("input", fp), 
            dtype=dtype_dict
        )
        for fp in os.listdir("input")
        if fp.endswith(suffix)
    ]).to_csv(
        suffix,
        index=None
    )