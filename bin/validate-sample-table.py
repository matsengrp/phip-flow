#!/usr/bin/env python

import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-s", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()

sample_table = pd.read_csv(args.s, header=0)
if "fastq_filepath" not in sample_table.columns:
    raise KeyError("Must provide 'fastq_filepath' column in the table")

if "sample_id" not in sample_table.columns:
    sample_table["sample_id"] = list(range(len(sample_table)))
    sample_table.set_index("sample_id", inplace=True)
else:
    sample_table.set_index("sample_id", inplace=True)

sample_table.to_csv(args.o, index=True, na_rep="NA")

