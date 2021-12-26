#!/usr/bin/env python
"""
convert peptide metadata to fasta format.
"""

import pandas as pd
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-pt", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()

def trim_index(sequence):
    return "".join([nt for nt in sequence if nt.isupper()])

fasta_fp = open(f"{sys.argv[1]}.fasta", "w")
with open(args.o, "w") as fasta_fp:
    peptide_table = pd.read_csv(args.pt, index_col=0, header=0)
    for index, row in peptide_table.iterrows():
        ref_sequence = trim_index(row["oligo"])
        fasta_fp.write(f">{index}\n{ref_sequence}\n")
