# PHIP-FLOW
A Nextflow pipeline for Common Phage Immuno-Precipitation Sequencing experiments.
See the [Documentation](https://matsengrp.github.io/phippery/introduction.html)
for more details and usage examples.

## Quickstart 

Install `Nextflow` by using the following command: 

    curl -s https://get.nextflow.io | bash 
    
Download the `Docker` Desktop, there exists several distributions packaged for
various linux flavors

    curl -fsSL https://get.docker.com -o get-docker.sh && sudo sh get-docker.sh

Launch the pipeline execution with the following command: 

    nextflow run matsengrp/phip-flow -r V1.12 -profile docker

Note: the ``-r VX.XX`` command runs the specified stable release version of the pipeline. 
For running the bleeding edge (not generally recommended) you may also specify ``-r main``.
You may also specify any of the 
[parameters](https://matsengrp.github.io/phippery/alignments-pipeline.html#parameters) 
for changing the input data and workflow behavior.

Note: the ``phippery`` [Dockerfile](https://github.com/matsengrp/phippery/blob/main/Dockerfile) 
contains all the required dependencies except those for ``EdgeR`` and ``BEER``, 
for which the maintainers of that package host their own public image.

[![Docker Repository on Quay](https://quay.io/repository/hdc-workflows/phippery/status "Docker Repository on Quay")](https://quay.io/repository/hdc-workflows/phippery)

