#!/usr/bin/env python
"""
This script should take in a dataset, sum the raw counts 
from all replicate sequences in the library, then proceed to
set the value for each replicate to that sum

Currently, this function only sets the raw counts, in place.
"""

import pandas as pd
import numpy as np
import phippery
from phippery.utils import *
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-ds", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()

ds = phippery.load(args.ds)

# find all value counts greater than 1,
pep_anno_table = get_annotation_table(ds, "peptide")

# Iterate over every group of peptides which share the same oligo sequence
for oligo_seq, pep_anno_table_oligo in pep_anno_table.groupby('oligo'):

    # Check to see if there are multiple peptides with the same oligo sequence
    if pep_anno_table_oligo.shape[0] == 1:

        # Don't make any changes for unique oligos
        continue

    # Otherwise, get the sum of the counts across all oligos
    idxs = pep_anno_table_oligo.index.values
    # rep_pep_sums = ds.counts.loc[idxs, :].sum(axis=0).values

    # Set the summed value for all peptides which share the same oligo sequence
    ds.counts.loc[idxs, :] = np.tile(
            ds.counts.loc[idxs, :].sum(axis=0).values,
            (pep_anno_table_oligo.shape[0], 1)
    )

phippery.dump(ds, args.o)
