# PhIP-Flow

A general PhIP-Seq short read alignment pipeline using [NextFlow](https://www.nextflow.io/) to produce a dataset containing the _raw alignment counts_ , _sample metadata table,_ and _peptide metadata table_ , all merged into an xarray [DataSet](http://xarray.pydata.org/en/stable/api.html#dataset). This dataset organization can subsequently be queried analyzed by using the [phippery](https://github.com/matsengrp/phippery) Python package.

For a usage template / example on simulated data (as well as a much more in-depth description of necessary configuration and input files) - see the [phip-flow-template](https://github.com/matsengrp/phip-flow-template).

The pipeline Directed Acyclic Graph is visualized below:
<p align="center">
  <img src="dag.png" width="375">
</p>
