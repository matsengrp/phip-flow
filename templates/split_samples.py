#!/usr/bin/env python3

from collections import defaultdict
import os
from typing import List
import pandas as pd
import logging


def setup_logging() -> logging.Logger:
    """Set up logging."""

    # Set up logging
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [split_samples] %(message)s'
    )
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Also write to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    return logger

logger = setup_logging()

# The user must specify a CSV containing the sample mapping
sample_mapping_fp = "!{params.dataset_prefix}_sample_annotation_table.csv.gz"
logger.info(f"Reading in sample mapping from: {sample_mapping_fp}")
assert os.path.exists(sample_mapping_fp)

# Read in the table
df = pd.read_csv(sample_mapping_fp, index_col=0)
logger.info(f"Sample mapping table has {df.shape[0]:,} rows and {df.shape[1]:,} columns")

# If the user specified a column used to group replicates
# from the same sample
sample_grouping_col = "!{params.sample_grouping_col}"
if len(sample_grouping_col) > 0:

    # Make sure that the column is present in the table
    msg = f"Column '{sample_grouping_col}' not found ({', '.join(df.columns.values)})"
    assert sample_grouping_col in df.columns.values, msg

    # Write out a file containing the unique list of sample names
    df.reindex(
        columns=[sample_grouping_col]
    ).drop_duplicates(
    ).to_csv(
        "sample_list",
        header=None,
        index=None
    )

# If no such grouping was found
else:

    # Just write out a list of each replicate
    with open("sample_list", "w") as handle:
        handle.write(
            "\n".join(list(map(str, df.index.values)))
        )
