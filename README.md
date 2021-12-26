# PHIP-FLOW
A Nextflow pipeline for Common Phage Immuno-Precipitation Sequencing experiments.
See the [Documentation](https://matsengrp.github.io/phippery/introduction.html)
for more details and usage examples.

[![nextflow]()]()
[![Build Status]()]()

## Quickstart 

Install `Nextflow` by using the following command: 

    curl -s https://get.nextflow.io | bash 
    
Download the `Docker` Desktop, there exists several distibutions packaged for
various linux flavors

    curl -fsSL https://get.docker.com -o get-docker.sh && sudo sh get-docker.sh

Launch the pipeline execution with the following command: 

    nextflow run matsengrp/phip-flow -profile docker

Note: the [Dockerfile](docker/Dockerfile) contains all the required dependencies. 
Add the `-profile docker` to enable the containerised execution to the 
example command line shown below. 

## Pipeline Description

## Input files

## Pipeline parameters
    
## Pipeline results

## Schematic Outline

## Requirements 

* [Nextflow](https://www.nextflow.io) 20.07.1 (or later)
* [Docker](https://www.docker.com/) 1.10 (or later) or [Singularity](http://singularity.lbl.gov) engine

required software components reported in the following section. See the included 
[Dockerfile](docker/Dockerfile) for the configuration details.

## Components 

PhIP-Flow uses the following software components and tools: 

<TODO>

## Acknowledgements

<TODO>
