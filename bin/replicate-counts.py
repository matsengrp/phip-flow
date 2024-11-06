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
from phippery.utils import load, dump, get_annotation_table
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-ds", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()


def replicate_oligo_counts(ds, peptide_oligo_feature="Oligo"):
    """This function should take in a dataset, sum the raw counts 
    from all replicate sequences in the library, then proceed to
    set the value for each replicate to that sum

    Currently, this function only sets the raw counts, in place.
    """

    # TODO remove this code
    # find all value counts greater than 1,
    #pep_anno_table = get_annotation_table(ds, "peptide")
    #oligo_vc = pep_anno_table["Oligo"].value_counts()

    ## for each oligo that is not unique in a library
    #for oligo, count in oligo_vc[oligo_vc > 1].items():
    #    replicate_idxs = pep_anno_table[
    #            pep_anno_table["Oligo"]==oligo
    #    ].index.values

    #    # sum the replicate values
    #    rep_pep_sums = ds.counts.loc[replicate_idxs, :].sum(axis=0).values

    #    # set the replicate counts equal to the sum of all
    #    ds.counts.loc[replicate_idxs, :] = np.tile(rep_pep_sums, (count, 1))

    # find all value counts greater than 1,
    pep_anno_table = get_annotation_table(ds, "peptide")

    # Iterate over every group of peptides which share the same oligo sequence
    for oligo_seq, pep_anno_table_oligo in pep_anno_table.groupby(peptide_oligo_feature):

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

ds = phippery.load(args.ds)
replicate_oligo_counts(ds, "oligo")
phippery.dump(ds, args.o)
