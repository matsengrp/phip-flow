#!/usr/bin/env python
"""
fit the zscore model.

Requirements:
1.  a sample table annotation column "control_status"
    with each sample either having the string factor level
    "beads_only" being the samples which we fit the model
    to, or "empirical" being the samples we predict on
    after the model is fit to each peptide.

2.  The xarray phip dataset passed in must have the
    "cpm" layer in the enrichment tables.
    and we expect the two types of sample groups
    were normalized using counts per million together
    as is the default when computing stats.

For a more complete description, 
please see the overview by Kevin Sung found at
https://matsengrp.github.io/phippery/
"""

from phippery.utils import * 
from phippery.modeling import zscore

import argparse
import warnings

parser = argparse.ArgumentParser()
parser.add_argument("-ds", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()

ds = load(args.ds)
beads_ds = ds_query(ds, "control_status == 'beads_only'")

zscore_ds = zscore(
    ds,
    beads_ds,
    data_table='cpm',
    min_Npeptides_per_bin=300,
    lower_quantile_limit=0.05,
    upper_quantile_limit=0.95,
    inplace=False,
    new_table_name='zscore'
)

dump(zscore_ds, args.o)
