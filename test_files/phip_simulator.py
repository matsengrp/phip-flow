"""
@Author: Jared Galloway
@date: 

phip-seq data set simulator.

The idea here being the simulation 
gives us an expected outcome from the pipeline,
this is purely for testing.
"""

import numpy as np
import os


class PhipSimulator(object):
    """
    simulates and holds onto metadata (is this a good place for @property?)
    so that reads files can be continuously added to.
    """

    def __init__(
        self,
    ):
        
        self.index_oligo_in_ref = {}
        self.sam_fq_fp = {}
        self.added_counts = None
        self.nt = ["A","T","C","G"]

    def generate_metadata(
        self,
        exp,
        ref,
        n_samples = 10,
        n_peptides = 10,
        fq_prefix = "sample-1",
        fq_pattern = "sample-*",
        virus="BANANA",
        protein="crazy_buffalo",
        tile_length = 93,
        adapt_3_length = 18,
        adapt_5_length = 16,
    ):
        """
        A general function to simulate random peptides, and create the sample
        metadata. Writes to f{lib}/sample_metadata.csv 
        and f{lib}/peptide_metadata.csv
        """

        # create peptide metadata by generating random nts
        ref_fp = open(f"{ref}/peptide_metadata.csv","w")
        adapter_3 = ''.join(np.random.choice(self.nt, adapt_3_length)).lower()
        adapter_5 = ''.join(np.random.choice(self.nt, adapt_5_length)).lower()
        ref_fp.write("ID,Virus,Protein,Loc,AA,Oligo,Peptide\n")
        for ID in range(n_samples):
            oligo = ''.join(np.random.choice(["A","T","C","G"], tile_length))
            self.index_oligo_in_ref[ID] = oligo
            oligo_w_adap = adapter_3 + oligo + adapter_5
            ref_fp.write(f"{ID},_,_,_,_,{oligo_w_adap},_\n")
        
        # create sample metadata
        sample_metadata_fp = open(f"samples/sample_metadata.csv","w")
        sample_metadata_fp.write("ID,fastq_pattern,experiment,reference,Sample_type,Notes\n")
        for sID in range(n_peptides):
            fq_name = f"{exp}/{fq_prefix}-{sID}.fastq"
            self.sam_fq_fp[sID] = open(fq_name,"w")
            sample_metadata_fp.write(f"{sID},{fq_pattern}-{sID}.fastq,{exp},{ref},_,_\n")
                

    def generate_reads(self, counts, n_mismatches=0, read_length=125):
        """ generate the fastq reads based off metadata """

        assert(self.index_oligo_in_ref != {})
        assert(self.sam_fq_fp != {})

        for pID in self.index_oligo_in_ref:
            for sID in self.sam_fq_fp:
                for read_idx in range(counts[pID][sID]):

                    # place the oligo from the reference somewhere in the
                    # middle of a read of length, read_length.
                    oligo = self.generate_mismatch(self.index_oligo_in_ref[pID], n_mismatches)
                    n_nt_filler = read_length - len(oligo)
                    filler = ''.join(np.random.choice(self.nt, n_nt_filler))
                    split = np.random.choice(range(n_nt_filler))
                    read = filler[:split] + oligo + filler[split:]
                    assert(len(read) == read_length)
                    qscore = "F" * read_length 
                    self.sam_fq_fp[sID].write(f"@\n{read}\n+\n{qscore}\n")


    def duplicate(self, source, destination):
        """ over write the oligo at position destination with source oligo 
        in the peptide metadata """
 
        assert source in self.index_oligo_in_ref.keys()
        assert destination in self.index_oligo_in_ref.keys()

        self.index_oligo_in_ref[destination] = self.index_oligo_in_ref[source]
        

    def generate_mismatch(self, sequence, num_mm):
        """
        Take a sequence of nucleotides and mutate (strictly to diff base)
        This should be used to emulate n mismatches during alignment.
        """

        nts = set(["A","T","C","G"])
        mismatch_sequence = ""
        mismatch_indices = np.random.choice(range(len(sequence)), num_mm, replace=False)
        for idx, nt in enumerate(sequence):
            if idx not in mismatch_indices:
                mismatch_sequence += nt
            else:
                mismatch_sequence += np.random.choice(list(nts - set(nt)))
        return mismatch_sequence
    

if __name__ == "__main__":
    
    np.random.seed(23)
    n_samples = 10
    n_peptides = 10

    ps = PhipSimulator()
    ps.generate_metadata(exp="lib_ones",ref="refa")

    # generate one "zero mismatch" read for each entry,
    # so the expected counts matrix will be 10 x 10 with 
    # a single hit for each sample peptide pair.
    no_mm = np.ones([n_samples,n_peptides]).astype(int)
    ps.generate_reads(no_mm, 0)

