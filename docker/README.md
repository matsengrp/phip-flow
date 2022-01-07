= Docker container

The pipline needs the following tools:

- http://www.htslib.org/[samtools]

A Docker container with all tools except the Genome Analysis Toolkit can be built from the `Dockerfile` present in this folder.

A prebuilt Docker container is available at the Docker hub as `cbcrg/ngs2017ws-nf`.

For the Genome Analysis Toolkit we use `matsengrp/phippery`.
