#!/usr/bin/env python
"""
"""

import pandas as pd
import numpy as np
import phippery
from phippery.utils import *
import argparse
import glob
import os
from functools import reduce
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-st", type=str)
parser.add_argument("-pt", type=str)
parser.add_argument("-cfp", type=str)
parser.add_argument("-sfp", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()


def _collect_sample_table(sample_table_filename: str):
    """Read and verify a sample table."""

    sample_table = pd.read_csv(sample_table_filename, sep=",", index_col=0, header=0)

    if sample_table.index.name != "sample_id":
        raise ValueError("The name of the index must be 'sample_id'")

    if sample_table.index.dtype != "int64":
        raise ValueError("The index values for sample_id must be inferred as integers")

    sample_table.sort_index(inplace=True)
    return sample_table


def _collect_peptide_table(peptide_table_filename: str):
    """Read and verify a peptide table."""

    peptide_table = pd.read_csv(peptide_table_filename, sep=",", index_col=0, header=0)

    if peptide_table.index.name != "peptide_id":
        raise ValueError

    if peptide_table.index.dtype != "int64":
        raise ValueError("The index values for peptide_id must be inferred as integers")

    peptide_table.sort_index(inplace=True)
    return peptide_table


def load_from_counts_tsv(
    sample_table, 
    peptide_table, 
    counts_file_pattern, 
    stats_file_pattern, 
):

    counts = [f for f in glob.glob(counts_file_pattern)]
    stats_files = [f for f in glob.glob(stats_file_pattern)]

    merged_counts = collect_counts(counts)
    peptide_table = _collect_peptide_table(peptide_table)
    sample_table = _collect_sample_table(sample_table)

    # For testing FIXME
    print(merged_counts)
    print(sample_table)
    
    # Add more summary metrics per sample
    sample_table = sample_table.assign(
        percent_mapped=sample_table["reads_mapped"] / sample_table["raw_total_sequences"] * 100.,
        percent_peptides_detected=(merged_counts > 0).mean() * 100.,
        percent_peptides_between_10_and_100=merged_counts.applymap(lambda v: v >= 10 & v <= 100).mean() * 100.,
    )
    
    # For testing FIXME
    print(sample_table)

    def num(s):
        try:
            return int(s)
        except ValueError:
            return float(s)

    if stats_files is not None:
        alignment_stats = defaultdict(list)
        for sample_alignment_stats in stats_files:
            fp = os.path.basename(sample_alignment_stats)
            sample_id = int(fp.strip().split(".")[0])
            alignment_stats["sample_id"].append(sample_id)
            for line in open(sample_alignment_stats, "r"):
                line = line.strip().split("\t")
                x = line[0]
                anno_name = "_".join(x.lower().split()).replace(":", "")
                alignment_stats[f"{anno_name}"].append(num(line[1]))

        stats_df = pd.DataFrame(alignment_stats).set_index("sample_id")
        
        sample_table = sample_table.merge(
                stats_df, 
                "outer", 
                left_index=True, 
                right_index=True
        )

    ds = stitch_dataset(
        counts=merged_counts, 
        peptide_table=peptide_table, 
        sample_table=sample_table,
    )

    #dump(ds, output)
    return ds


ds = load_from_counts_tsv(
    args.st,
    args.pt,
    args.cfp,
    args.sfp,
)

dump(ds, args.o)

#def merge_count_data(counts):
#    """
#    This function takes in a list of paths which
#    contains the counts for each peptide alignment
#    for each sample. These files should contain
#    no header.
#
#    :param: counts <str> - a list of paths leading
#    to raw peptide enrichment counts for each sample
#    """
#
#    load = lambda path, sample: pd.read_csv(  # noqa
#        path, index_col=0, sep="\t", names=["peptide_id", sample]
#    )
#
#    sample_dataframes = [
#        load(path, int(os.path.basename(path).split(".")[0])) for path in counts
#    ]
#
#    merged_counts_df = reduce(
#        lambda l, r: pd.merge(l, r, how="outer", left_index=True, right_index=True),
#        sample_dataframes,
#    ).fillna(0)
#
#    merged_counts_df.columns = merged_counts_df.columns.astype(int)
#    merged_counts_df.index = merged_counts_df.index.astype(int)
#    merged_counts_df.sort_index(inplace=True)
#    merged_counts_df.sort_index(axis=1, inplace=True)
#
#    return merged_counts_df


#def load_from_counts_tsv(
#    sample_table, 
#    peptide_table, 
#    counts_file_pattern, 
#    stats_file_pattern, 
#    output
#):
#    """
#    Collect sample and peptide metadata tables along with a
#    two-column tsv file for each sample, 
#    and produce a properly formatted xarray dataset.
#    """
#
#    counts = [f for f in glob.glob(counts_file_pattern)]
#    stats_files = [f for f in glob.glob(stats_file_pattern)]
#
#    merged_counts = collect_counts(counts)
#    peptide_table = collect_peptide_table(peptide_table)
#    sample_table = collect_sample_table(sample_table)
#
#    def num(s):
#        try:
#            return int(s)
#        except ValueError:
#            return float(s)
#
#    if stats_files is not None:
#        alignment_stats = defaultdict(list)
#        for sample_alignment_stats in stats_files:
#            fp = os.path.basename(sample_alignment_stats)
#            sample_id = int(fp.strip().split(".")[0])
#            alignment_stats["sample_id"].append(sample_id)
#            for line in open(sample_alignment_stats, "r"):
#                line = line.strip().split("\t")
#                x = line[0]
#                anno_name = "_".join(x.lower().split()).replace(":", "")
#                alignment_stats[f"{anno_name}"].append(num(line[1]))
#
#        stats_df = pd.DataFrame(alignment_stats).set_index("sample_id")
#        
#        sample_table = sample_table.merge(
#                stats_df, 
#                "outer", 
#                left_index=True, 
#                right_index=True
#        )
#
#    ds = stitch_dataset(
#        counts=merged_counts, 
#        peptide_table=peptide_table, 
#        sample_table=sample_table,
#    )
#
#    dump(ds, output)

