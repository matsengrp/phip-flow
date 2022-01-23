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

import phippery
import phippery.modeling as modeling
import phippery.phipdata as phipdata
import phippery.utils as utils
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-ds", type=str)
parser.add_argument("-o", type=str)
args = parser.parse_args()

ds = phippery.load(args.ds)

# grab the relevant mock ip samples id's
beads_ids = utils.sample_id_coordinate_from_query(
        ds, ["control_status == 'beads_only'"]
)

beads_ds = ds.loc[dict(sample_id=beads_ids)]

params, fit_ds = modeling.neg_binom_model(
    ds,
    beads_ds,
    nb_p=2,
    trim_percentile=100.,
    outlier_reject_scale=10.,
    data_table="size_factors",
    inplace=False,
    new_table_name="neg_binom_mlxp",
)

phippery.dump(fit_ds, args.out)
