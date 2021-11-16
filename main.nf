
// TODO - remove automatic publish_dir from processes generate counts, leave that to io
// TODO - remove workflow from enrichment and just write some generic processes
// TODO write a batch worflow, which parallelizes the given workflow around sample 
// (or peptide? groups)
// TODO add a check workflow that makes sure the columns for requested workflows exist
// could make a replacement config file - meh later
// TODO add a test of example data
// TODO separate alignment into template which can be easily manipulated if desired
// TODO clean up containers
// TODO: can we split these things up? like a config for infrastructure?
// TODO: move the phippery workflow to processes and out of templates - keep it simple

nextflow.enable.dsl=2
include { ALIGN_COUNTS as generate_alignment_counts } from './workflows/alignment-counts-workflow.nf'
include { cpm_fold_enrichment as compute_cpm_fold_enrichment } from './workflows/enrichment-workflow.nf'
include { to_tall as write_ds_to_tall_csv } from './workflows/io.nf'
include { to_wide as write_ds_to_wide_csv } from './workflows/io.nf'

workflow {
    // GENERATE ALIGNMENT COUNTS
    ds = generate_alignment_counts()

    // ANALYSIS
    if( params.compute_enrichment )
        ds = compute_cpm_fold_enrichment(ds)
   
    // IO
    if( params.output_tall ) 
        write_ds_to_tall_csv(ds)
    if( params.output_wide )
        write_ds_to_wide_csv(ds)
}
