# PhIP-Flow

A general PhIP-Seq short read alignment pipeline using [NextFlow](https://www.nextflow.io/) to produce a dataset containing the _raw alignment counts_ , _sample metadata table,_ and _peptide metadata table_ , all merged into an xarray [DataSet](http://xarray.pydata.org/en/stable/api.html#dataset). This dataset organization can subsequently be queried analyzed by using the [phippery](https://github.com/matsengrp/phippery) Python package.

This repository exists primarily to host the bleading edge pipeline script, and should not need to be cloned unless
one would like to modify the pipeline. If looking to simply run with the appropriate 
confguration scripts, 
one can simply use the Nextflow's build in 
git aware infrastructure
```
$ nextflow run matsengrp/phip-flow/PhIP-Flow.nf -C foo-phip.config -o bar.phip
```

For a usage template / example on simulated data (as well as a much more in-depth description of necessary configuration and input files) - see the [phip-flow-template](https://github.com/matsengrp/phip-flow-template).

The pipeline Directed Acyclic Graph is visualized below:
<p align="center">
  <img src="dag.png" width="375">
</p>

For example configuration files and more helpful materials, the 
[Documentation](https://matsengrp.github.io/phippery/Nextflow.html)
will walk you through examples 
that can be found in the 
[Template repository](https://github.com/matsengrp/phip-flow-template)
