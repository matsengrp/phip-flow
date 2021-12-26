#!/usr/bin/env python

import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-p", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()

peptide_table = pd.read_csv(args.p, header=0)
if "oligo" not in peptide_table.columns:
    raise KeyError("Must provide 'oligo' column in the table")

if "peptide_id" not in peptide_table.columns:
    peptide_table["peptide_id"] = list(range(len(peptide_table)))
    peptide_table.set_index("peptide_id", inplace=True)
else:
    peptide_table.set_index("peptide_id", inplace=True)
    
peptide_table.to_csv(args.o, index=True, na_rep="NA")
