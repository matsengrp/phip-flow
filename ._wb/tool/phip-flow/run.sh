#!/bin/bash

set -meuo pipefail

date
echo
echo "Running workflow from ${PWD}"
echo

# Run the workflow
echo Starting workflow
nextflow \
    run \
    "${TOOL_REPO}" \
    --results "${PWD}" \
    -params-file ._wb/tool/params.json \
    -resume \
    -profile "${PROFILE}" &

# Get the process ID
PID="$!"

# Make a task which can kill this process
if [ ! -d ._wb/bin ]; then mkdir ._wb/bin; fi
echo """
#!/bin/bash

echo \"\$(date) Sending a kill signal to the workflow\"

kill ${PID}
""" > ._wb/bin/stop
chmod +x ._wb/bin/stop

# Bring the command back to the foreground
fg %1

echo
date
echo Done
