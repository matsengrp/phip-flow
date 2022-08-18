#!/usr/bin/env python
"""
Fit a negative binomial model to beads_only samples in
a dataset. Then predict on empirical IP samples to
get mlxp values of significance.

Requirements:
1.  a sample table annotation column "control_status"
    with each sample either having the string factor level
    "beads_only" being the samples which we fit the model
    to, or "empirical" being the samples we predict on
    after the model is fit to each peptide.

2.  The xarray phip dataset passed in must have the
    "size_factors" layer in the enrichment tables.
    and we expect the two types of sample groups
    were normalized using size factors together
    as is the default when computing stats.

For a more complete description, 
please see the overview by Kevin Sung found at
https://matsengrp.github.io/phippery/bkg-model.html 
"""

from phippery.modeling import neg_binom_model
from phippery.utils import *

import argparse
import warnings

parser = argparse.ArgumentParser()
parser.add_argument("-ds", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()

ds = load(args.ds)
beads_ds = ds_query(ds, "control_status == 'beads_only'")

if len(beads_ds.sample_id) <= 50:
    warnings.warn(
    "With less that 50 beads_only samples, "
    "it's probably that many peptide distributions fits "
    "will not converge, esspecially if coverage is low. "
    "See https://matsengrp.github.io/phippery/bkg-model.html for more."
    )

params, fit_ds = neg_binom_model(
    ds,
    beads_ds,
    nb_p=2,
    trim_percentile=100.,
    outlier_reject_scale=10.,
    data_table="size_factors",
    inplace=False,
    new_table_name="neg_binom_mlxp",
)

dump(fit_ds, args.o)
