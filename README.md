# PhIP-Flow

An alignment pipeline to produce an organized and coherent phip seq dataset, which can subsequently be used by the [phippery](https://github.com/matsengrp/phippery) Python package for query, analysis, and more.

This repository exists primarily to host the bleading edge pipeline script, and should not need to be cloned unless
one would like to modify the pipeline. If looking to simply run with the appropriate 
confguration scripts, 
one can simply use the Nextflow's build in 
git aware infrastructure
```
$ nextflow run matsengrp/phip-flow/PhIP-Flow.nf -C foo-phip.config -o bar.phip
```

The pipeline Directed Acyclic Graph is visualized below
<p align="center">
  <img src="dag.png" width="375">
</p>

For example configuration files and more helpful materials, the 
[Documentation](https://matsengrp.github.io/phippery/Nextflow.html)
will walk you through examples 
that can be found in the 
[Template repository](https://github.com/matsengrp/phip-flow-template)
