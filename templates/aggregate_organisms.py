#!/usr/bin/env python3

import os
from typing import List
import pandas as pd
import logging
from scipy.stats import gmean

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
# set of epitopes which do not overlap with any other epitope by more than 7aa.
# To identify overlaps, use an exact alignment approach using k-mers. Note that this
# will count two peptides as overlapping if they share any 7aa sequence, without
# performing global alignment.

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

        # The user must specify the maximum overlap
        self.max_overlap = int("!{params.max_overlap}")
        self.logger.info(f"Maximum overlap: {self.max_overlap}")

        # The user must specify the minimum z-score threshold
        self.zscore_threshold = float("!{params.zscore_threshold}")
        self.logger.info(f"Z-score threshold: {self.zscore_threshold}")

        # Read in the z-scores
        zscores_fp = "!{params.dataset_prefix}_zscore.csv.gz"
        self.logger.info(f"Reading in z-scores from: {zscores_fp}")
        assert os.path.exists(zscores_fp)
        self.zscores = pd.read_csv(zscores_fp, index_col=0)

        # Read in the edgeR hits (if present)
        edgeR_hits_fp = "!{params.dataset_prefix}_edgeR_hits.csv.gz"
        if os.path.exists(edgeR_hits_fp):
            self.logger.info(f"Reading in edgeR hits from: {edgeR_hits_fp}")
            self.edgeR_hits = pd.read_csv(
                edgeR_hits_fp,
                index_col=0,
                true_values=["TRUE", "True", "true"],
                false_values=["FALSE", "False", "false"],
                na_values=["NA", "N/A", "Na", "na", "n/a"]
            )
            self.has_edgeR_hits = True
        else:
            self.has_edgeR_hits = False

        # Group the replicates by sample
        self.logger.info("Grouping replicates by sample")
        self.sample_table = self.group_replicates()

        # Apply the max_overlap filter
        # (setting the column 'passes_filter' to True if the peptide passes)
        self.sample_table = self.apply_max_overlap_filter()

        # Save to CSV
        self.sample_table.to_csv("!{sample_id}.peptide.ebs.csv.gz", index=None)

        # Group the peptides by organism
        self.logger.info("Grouping peptides by organism")
        self.organism_table = self.group_organisms()

        # Save to CSV
        self.logger.info("Writing organism-level outputs to CSV")
        self.organism_table.to_csv("!{sample_id}.organism.summary.csv.gz", index=None)

        self.logger.info("Done")

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
        sample_mapping_fp = "!{params.dataset_prefix}_sample_annotation_table.csv.gz"
        self.logger.info(f"Reading in sample mapping from: {sample_mapping_fp}")
        assert os.path.exists(sample_mapping_fp)

        # Read in the table
        df = pd.read_csv(sample_mapping_fp, index_col=0)
        self.logger.info(f"Sample mapping table has {df.shape[0]:,} rows and {df.shape[1]:,} columns")

        # If the user specified a column used to group replicates
        # from the same sample
        sample_grouping_col = "!{params.sample_grouping_col}"
        if len(sample_grouping_col) > 0:

            # Make sure that the column is present in the table
            msg = f"Column '{sample_grouping_col}' not found ({', '.join(df.columns.values)})"
            assert sample_grouping_col in df.columns.values, msg

            # Return the column mapping of replicates to samples
            return df[sample_grouping_col]

        # Otherwise, if no grouping was specified
        else:

            # Just treat each sample the same
            return {
                int(replicate_id): str(replicate_id)
                for replicate_id in df.index.values
            }

    def read_peptide_mapping(self) -> pd.DataFrame:
        """Read the table mapping peptides (by ID) to organism, protein, and start position ('pos')."""

        peptide_mapping_fp = "!{params.dataset_prefix}_peptide_annotation_table.csv.gz"
        self.logger.info(f"Reading in peptide mappings from: {peptide_mapping_fp}")
        assert os.path.exists(peptide_mapping_fp)

        # Read in the table
        df = pd.read_csv(peptide_mapping_fp, index_col=0)
        self.logger.info(f"Peptide mapping table has {df.shape[0]:,} rows and {df.shape[1]:,} columns")

        # Map the user-provided names to controlled values
        mapping = {
            # The user must specify the column used to group peptides by organism
            "!{params.peptide_org_col}": "organism",
            # And by the protein sequence (which corresponds to the public epitope sequences)
            "!{params.peptide_seq_col}": "seq"
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

        # Add the length of the peptide
        df = df.assign(
            peptide_length=lambda d: d["seq"].apply(len)
        )

        return df

    def read_public_epitopes(self) -> List[str]:
        """Read the list of public epitopes provided."""

        # Table of public epitopes
        df = pd.read_csv("!{public_epitopes_csv}")
        self.logger.info(f"Public epitope table has {df.shape[0]:,} rows")

        # The user must specify the column which contains the public epitopes
        public_epitopes_col = "peptide_translate"

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

        # Get the replicates which should be combined for this sample
        replicates = [
            rep_i
            for rep_i in self.zscores.columns.values
            if self.sample_mapping.get(int(rep_i)) == '!{sample_id}'
        ]

        self.logger.info(f"Filtering down to the {len(replicates):,} replicates for sample '!{sample_id}'")
        assert len(replicates) > 0

        # Take a slice of the table
        df = self.zscores.reindex(columns=replicates)
        # If we have edgeR data, filter it to the same replicates
        if self.has_edgeR_hits:
            self.edgeR_hits = self.edgeR_hits.reindex(columns=replicates)
            
        # Add summary metrics
        df = df.assign(
            n_replicates=len(replicates),
            EBS=df.mean(axis=1),
            hit=df.apply(self.classify_hit, axis=1),
            edgeR_hit=(
                self.edgeR_hits.apply(self.classify_edgeR_hit, axis=1)
                if self.has_edgeR_hits
                else None
            ),
            sample='!{sample_id}'
        ).reset_index(
        ).rename(
            columns=dict(index="peptide")
        ).drop(
            columns=replicates + (
                ["edgeR_hit"] if not self.has_edgeR_hits else []
            )
        )

        # Mark whether each peptide is public
        df = df.assign(
            public=df["peptide"].apply(int).apply(
                lambda i: self.peptide_mapping["public"][i]
            )
        )

        return df

    def classify_hit(self, r):
        """Determine whether a peptide is a hit, or discordant."""

        # Get the vector of whether each replicate
        # is above the z-score threshold
        hit_vec = r > self.zscore_threshold

        # Determine the hit type
        if hit_vec.all():
            return "TRUE"
        elif not hit_vec.any():
            return "FALSE"
        else:
            return "DISCORDANT"

    def classify_edgeR_hit(self, r: pd.Series) -> str:
        """
        Determine whether a peptide is a hit, or discordant - 
        based on edgeR hits.
        """

        # Drop NA values before classification
        r = r.dropna()
        if len(r) == 0:  # If all values were NA
            return "NA"

        # Determine the hit type
        if r.all():
            return "TRUE"
        elif not r.any():
            return "FALSE"
        else:
            return "DISCORDANT"

    def apply_max_overlap_filter(self) -> pd.DataFrame:
        """Apply the max_overlap filter to each sample/organism."""

        # Analyze each sample/organism independently
        df = pd.concat([
            self.apply_max_overlap_filter_sub(d)
            for _, d in self.sample_table.assign(
                organism=lambda d: d["peptide"].apply(
                    self.peptide_mapping["organism"].get
                )
            ).groupby(
                ["sample", "organism"]
            )
        ])

        return df

    def apply_max_overlap_filter_sub(
        self,
        df: pd.DataFrame
    ) -> pd.DataFrame:

        # Add the sequence information for each peptide
        df = df.assign(
            seq=df["peptide"].apply(
                self.peptide_mapping["seq"].get
            ).apply(
                lambda s: s.rstrip("*")
            )
        )

        # Sort by EBS (descending)
        df = df.sort_values(by="EBS", ascending=False)

        # Keep track of the peptide kmers which have been observed so far
        kmers_seen = set()

        # Make a list of the indices pass the filter
        passes_filter = list()

        # Go down the list, starting with the tightest binders
        for _, r in df.iterrows():

            # Get the kmers by this peptide
            row_kmers = set([
                r["seq"][n:(n + self.max_overlap)]
                for n in range(len(r["seq"]) - self.max_overlap)
            ])

            # If none of those kmers have been seen before,
            # it passes the filter
            passes_filter.append(len(row_kmers & kmers_seen) == 0)

            # If it passes
            if passes_filter[-1]:

                # Add the covered positions
                kmers_seen |= row_kmers

        # Add a column to the table indicating
        # whether the peptide passes the filter
        df = df.assign(
            passes_filter=passes_filter
        )

        # Drop the sequence column
        return (
            df
            .drop(columns=["seq"])
            .sort_index()
        )

    def group_organisms(self) -> pd.DataFrame:
        """Group together the results by organism."""

        # Analyze each organism independently
        df = pd.concat([
            self.group_sample_organisms(d, sample, organism)
            for (sample, organism), d in self.sample_table.assign(
                organism=lambda d: d["peptide"].apply(
                    self.peptide_mapping["organism"].get
                )
            ).groupby(
                ["sample", "organism"]
            )
        ]).fillna(
            0
        )

        return df

    def group_sample_organisms(
        self,
        df: pd.DataFrame,
        sample: str,
        organism: str
    ) -> pd.DataFrame:

        """Analyze the data for a single sample, single organism."""

        # For this summary, drop peptides which don't pass the filter
        df = df.query("passes_filter")

        # Return the number of hits, etc. for all and just public epitopes
        dat = pd.DataFrame([{
            "sample": sample,
            "organism": organism,
            **{
                k: v
                for label, d in [
                    ("all", df),
                    ("public", df.query("public")),
                    ("hits", df.query("hit == 'TRUE'")),
                ]
                if d.shape[0] > 0
                for k, v in [
                    (f"n_hits_{label}", (d["hit"] == "TRUE").sum()),
                    (f"n_discordant_{label}", (d["hit"] == "DISCORDANT").sum()),
                    (f"max_ebs_{label}", d["EBS"].max()),
                    (f"mean_ebs_{label}", d["EBS"].mean()),
                    (f"gmean_ebs_{label}", gmean(d["EBS"]))
                ] + (
                    [
                        (f"n_edgeR_hits_{label}", (d["edgeR_hit"] == "TRUE").sum()),
                        (f"n_edgeR_discordant_{label}", (d["edgeR_hit"] == "DISCORDANT").sum()),
                    ]
                    if self.has_edgeR_hits
                    else []
                )
                if k not in [
                    "n_hits_hits",
                    "n_discordant_hits",
                    "gmean_ebs_all",
                    "gmean_ebs_public"
                ]
            }
        }])

        return dat


AggregatePhIP()
