
nextflow.enable.dsl=2
include { ALIGN_COUNTS as generate_alignment_counts } from './workflows/alignment-counts-workflow.nf'
include { cpm_fold_enrichment as compute_cpm_fold_enrichment } from './workflows/enrichment-workflow.nf'
include { to_tall as write_ds_to_tall_csv } from './workflows/io.nf'
include { to_wide as write_ds_to_wide_csv } from './workflows/io.nf'

workflow {
    // GENERATE ALIGNMENT COUNTS

    // Input defined by the nextflow.config file
    // If you rename or make your own config file, 
    // be sure to specificy -C your.config when running
    // this workflow.

    // output 'ds' is the pickle dump'd binary xarray
    // dataset containing the raw counts
    ds = generate_alignment_counts()

    // ANALYSIS
    // the config file specifies which analysis you would like run.
    // each one of these subflows takes the binary xarray file,
    // and returns it modified ("layered") with the analysis
    if( params.compute_enrichment )
        ds = compute_cpm_fold_enrichment(ds)
   
    // OUTPUT
    // Write the analysis to disk in any of the formats you desire.
    if( params.output_tall ) 
        write_ds_to_tall_csv(ds)

    if( params.output_wide )
        write_ds_to_wide_csv(ds)
}
