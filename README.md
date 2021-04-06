# PhIP-Flow

A general PhIP-Seq short read alignment pipeline to produce a dataset conatining the raw alignment counts for each sample merged into an xarray [DataSet](). This dataset organization can subsequently be easily analyzed by using the the [phippery](https://github.com/matsengrp/phippery) Python package.

The pipeline Directed Acyclic Graph is visualized below
<p align="center">
  <img src="dag.png" width="375">
</p>

For a config template / example on simulated data - as well as a much more in-depth description of necessary configurations and input files - see the [phip-flow-template](https://github.com/matsengrp/phip-flow-template)
