#!/bin/bash

# Note: Setting the -m flag is essential for the fg command
set -eumo pipefail

echo "Setting up nextflow.config"

echo """
workDir = '${WORK_DIR}'
singularity.cacheDir = '${CACHE_DIR}'
docker.enabled = true
report.enabled = true
trace.enabled = true
""" > nextflow.config

cat nextflow.config
echo

# Disable ANSI logging
export NXF_ANSI_LOG=false

# Print the Nextflow version being used
echo "Nextflow Version: ${NXF_VER}"
echo

# Execute the tool in the local environment
echo "Starting tool"
echo

# Start the tool
/bin/bash ._wb/helpers/run_tool
