#!/usr/bin/env python3

from collections import defaultdict
import os
from typing import List
import pandas as pd
import logging

# APPROACH

# 1. Start with the Z-scores calculated per epitope across every sample replicate

# 2. Combine the epitope-level data collected for multiple
# replicates of the same sample, saving a table with the
# following information per epitope, per sample:
#   Mean Z-score, also referred to as the Epitope Binding Score (EBS)
#   Hit:
#       marked as TRUE if both replicates are above the threshold Z-score
#       marked as FALSE if both replicates are below the threshold Z-score
#       marked as DISCORDANT if some but not all replicates are above the threshold Z-score
#   Public: marked as TRUE if the epitope was included in the input list of public epitopes
 
# 3. To combine the virus-level data for each sample, only keep the highest-scoring
# set of epitopes which do not overlap with any other epitope by more than 7aa

# 4. Finally, save a table with the following information per virus, per sample:
#   Number of all epitope hits
#   Number of public epitope hits
#   Number of all discordant epitopes
#   Number of public discordant epitopes
#   Max EBS across all epitopes
#   Max EBS across public epitopes
#   Mean EBS across all epitopes
#   Mean EBS across public epitopes
 

# INPUTS

# Read in the input data from an expected location
# The placement of the appropriate files in these locations
# is expected to be performed by the Nextflow wrapper code
# using the configuration of the appropriate module within
# the phip-flow workflow

class AggregatePhIP:

    def __init__(self):

        # Set up logging
        self.logger = self.setup_logging()

        # Mapping replicates to samples
        self.sample_mapping = self.read_sample_mapping()

        # List of public epitopes
        self.public_epitopes = self.read_public_epitopes()

        # Mapping peptides to organisms
        self.peptide_mapping = self.read_peptide_mapping()

        # The user must specify the length of each peptide
        self.peptide_length = int("${params.peptide_length}")
        self.logger.info(f"Peptide length: {self.peptide_length}")

        # The user must specify the maximum overlap
        self.max_overlap = int("${params.max_overlap}")
        self.logger.info(f"Maximum overlap: {self.max_overlap}")

        # The user must specify the minimum z-score threshold
        self.zscore_threshold = float("${params.zscore_threshold}")
        self.logger.info(f"Z-score threshold: {self.zscore_threshold}")

        # Read in the z-scores
        zscores_fp = "${params.dataset_prefix}_zscore.csv"
        self.logger.info(f"Reading in z-scores from: {zscores_fp}")
        assert os.path.exists(zscores_fp)
        self.zscores = pd.read_csv(zscores_fp, index_col=0)

        # Group the replicates by sample
        self.sample_table = self.group_replicates()

        # Save to CSV
        self.sample_table.to_csv("peptide.ebs.csv.gz", index=None)

        # Group the peptides by organism
        self.organism_table = self.group_organisms()

        # Save to CSV
        self.organism_table.to_csv("organism.summary.csv.gz", index=None)

    def setup_logging(self) -> logging.Logger:
        """Set up logging."""

        # Set up logging
        logFormatter = logging.Formatter(
            '%(asctime)s %(levelname)-8s [aggregate_organisms] %(message)s'
        )
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)

        # Also write to STDOUT
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        logger.addHandler(consoleHandler)

        return logger

    def read_sample_mapping(self) -> pd.Series:
        """Read a mapping of replicates to samples."""

        # The user must specify a CSV containing the sample mapping
        sample_mapping_fp = "${params.dataset_prefix}_sample_annotation_table.csv"
        self.logger.info(f"Reading in sample mapping from: {sample_mapping_fp}")
        assert os.path.exists(sample_mapping_fp)

        # Read in the table
        df = pd.read_csv(sample_mapping_fp, index_col=0)
        self.logger.info(f"Sample mapping table has {df.shape[0]:,} rows and {df.shape[1]:,} columns")

        # The user must specify the column used to group replicates
        # from the same sample
        sample_grouping_col = "${params.sample_grouping_col}"

        msg = f"Column '{sample_grouping_col}' not found ({', '.join(df.columns.values)})"
        assert sample_grouping_col in df.columns.values, msg

        # Return the column mapping of replicates to samples
        return df[sample_grouping_col]

    def read_peptide_mapping(self) -> pd.DataFrame:
        """Read the table mapping peptides (by ID) to organism, protein, and start position ('pos')."""

        peptide_mapping_fp = "${params.dataset_prefix}_peptide_annotation_table.csv"
        self.logger.info(f"Reading in peptide mappings from: {peptide_mapping_fp}")
        assert os.path.exists(peptide_mapping_fp)

        # Read in the table
        df = pd.read_csv(peptide_mapping_fp, index_col=0)
        self.logger.info(f"Peptide mapping table has {df.shape[0]:,} rows and {df.shape[1]:,} columns")

        # Map the user-provided names to controlled values
        mapping = {
            # The user must specify the column used to group peptides by organism
            "${params.peptide_org_col}": "organism",
            # By protein name
            "${params.peptide_prot_col}": "protein",
            # By the starting position of the peptide within the protein
            "${params.peptide_pos_col}": "position",
            # And by the protein sequence (which corresponds to the public epitope sequences)
            "${params.peptide_seq_col}": "seq"
        }

        # For each of the user-provided columns
        for cname in mapping.keys():

            # Make sure that it is in the table
            msg = f"Column '{cname}' not found ({', '.join(df.columns.values)})"
            assert cname in df.columns.values, msg

        # Change the names
        df = df.rename(columns=mapping)

        # Only return those columns
        df = df.reindex(
            columns=list(mapping.values())
        )

        # Assign the column `public` True if the protein sequence is in the public epitope list
        # Make sure to strip everything after the "*"
        df = df.assign(
            public=df["seq"].apply(
                lambda s: s.split("*")[0]
            ).isin(self.public_epitopes)
        )

        self.logger.info(f"Public Epitopes: {df['public'].sum():,} / {df.shape[0]:,}")

        # Drop the sequence
        df = df.drop(
            columns=["seq"]
        )

        return df

    def read_public_epitopes(self) -> List[str]:
        """Read the list of public epitopes provided."""

        # Table of public epitopes
        df = pd.read_csv("${public_epitopes_csv}")
        self.logger.info(f"Public epitope table has {df.shape[0]:,} rows")

        # The user must specify the column which contains the public epitopes
        public_epitopes_col = "peptide_translate"                             # FIXME

        msg = f"Column not found: {public_epitopes_col} in ({', '.join(df.columns.values)})"
        assert public_epitopes_col in df.columns.values, msg

        # Strip everything after the "*"
        return df[
            public_epitopes_col
        ].apply(
            lambda s: s.split("*")[0]
        ).tolist()

    def group_replicates(self) -> pd.DataFrame:
        """Group together the replicates of the same sample."""

        # Make a long table and add the sample designation
        df = self.zscores.reset_index(
        ).melt(
            id_vars=["index"]
        ).rename(
            columns=dict(
                index="peptide",
                variable="replicate",
                value="zscore"
            )
        ).assign(
            sample=lambda d: d["replicate"].apply(int).apply(self.sample_mapping.get)
        ) #.query("replicate < 5") # FIXME

        # Make sure that we have sample labels for everything
        assert not df["sample"].isnull().any(), "Could not find sample labels for all replicates"

        # For each peptide, for each sample, calculate the summary metrics
        df = df.groupby(
            ["sample", "peptide"]
        ).apply(
            self.group_peptides_replicates
        ).reset_index(
        ).drop(
            columns=["level_2"]
        )
        
        # Mark whether each peptide is public
        df = df.assign(
            public=df["peptide"].apply(int).apply(
                lambda i: self.peptide_mapping["public"][i]
            )
        )

        return df

    def group_peptides_replicates(self, df) -> pd.DataFrame:

        # Get the vector of whether each replicate is above the z-score threshold
        hit_vec = df["zscore"] > self.zscore_threshold

        # Determine the hit type
        if hit_vec.all():
            hit = "TRUE"
        elif not hit_vec.any():
            hit = "FALSE"
        else:
            hit = "DISCORDANT"

        return pd.DataFrame(
            [
                dict(
                    n_replicates=df.shape[0],
                    EBS=df["zscore"].mean(),
                    hit=hit
                )
            ]
        )

    def group_organisms(self) -> pd.DataFrame:
        """Group together the results by organism."""

        # Analyze each organism independently
        df = self.sample_table.assign(
            organism=lambda d: d["peptide"].apply(
                self.peptide_mapping["organism"].get
            )
        ).groupby(
            ["sample", "organism"]
        ).apply(
            self.group_sample_organisms
        ).reset_index(
        ).fillna(
            0
        )

        return df

    def group_sample_organisms(self, df:pd.DataFrame) -> pd.DataFrame:
        """Analyze the data for a single sample, single organism."""

        # Add the protein and position labels
        df = df.assign(
            protein=df["peptide"].apply(
                self.peptide_mapping["protein"].get,
            ),
            position=df["peptide"].apply(
                self.peptide_mapping["position"].get,
            )
        )

        # Sort by EBS (descending)
        df = df.sort_values(by="EBS", ascending=False)

        # Keep track of which positions have been covered
        covered_positions = defaultdict(set)

        # Make a list of the indices which will be dropped
        to_drop = list()

        # Go down the list, starting with the tightest binders
        for i, r in df.iterrows():

            # Get the positions covered by this peptide
            row_pos = set(range(r["position"], r["position"] + self.peptide_length))

            # Get the number of overlapping positions
            n_overlap = len(covered_positions[r["protein"]] & row_pos)

            # If the maximum overlap threshold is exceeded
            if n_overlap >= self.max_overlap:

                # Drop the row
                to_drop.append(i)

            # If not
            else:

                # Add the covered positions
                covered_positions[r["protein"]] |= row_pos

        df = df.drop(index=to_drop)

        # Return the number of hits, etc. for all and just public epitopes
        return pd.DataFrame([{
            k: v
            for label, d in [
                ("all", df),
                ("public", df.query("public")),
            ]
            if d.shape[0] > 0
            for k, v in [
                (f"n_hits_{label}", (d["hit"] == "TRUE").sum()),
                (f"n_discordant_{label}", (d["hit"] == "DISCORDANT").sum()),
                (f"max_ebs_{label}", d["EBS"].max()),
                (f"mean_ebs_{label}", d["EBS"].mean())
            ]
        }])


AggregatePhIP()