#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// This pipeline will compare predicted resistance genes with measured antibiotic resistance
// to identify genes that may be associated with carbapenem resistance. The outcome of this
// pipeline is a table that may be used to update the list of carbapenem associated genes in
// the typing pipeline.

// Container versions
r_container = "stitam/r-aci:0.13"

// Database versions
bvbrc_db_ver = "v230216"

// Output parameters
params.resdir = "${launchDir}/carbapenem"

// Processes

process calibrate_carbapenem_prediction {
    container "$r_container"
    storeDir "${params.resdir}/calibrate_carbapenem_prediction"

    input:
    path aci_with_qc
    path amr_bvbrc

    output:
    path "aci_wide.rds"
    path "amr_wide.rds"
    path "determinants.rds"
    path "determinants.tsv"
    path "multiple_determinants.rds"
    path "ISAba1.rds"
    path "precision_recall.pdf"
    path "roc.pdf"

    script:
    """
    Rscript $projectDir/bin/calibrate_carbapenem_prediction.R \
    --project_dir $projectDir \
    --file $aci_with_qc \
    --resgene_groups $projectDir/data/resgene_groups.csv \
    --amr_db $amr_bvbrc
    """
}

process run_qc_checks {
    container "$r_container"
    storeDir "${params.resdir}/run_qc_checks"

    output:
    path "aci_with_qc.rds", emit: aci_with_qc
    path "aci_with_qc.tsv"
    path "aci_qc_failed.rds"
    path "aci_qc_failed.tsv"
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/run_qc_checks.R \
    --project_dir $projectDir \
    --file $params.assembly_summary
    """
}

process tidy_bvbrc {
    container "$r_container"
    storeDir "${params.resdir}/tidy_bvbrc"

    output:
    path "amr_bvbrc.rds"

    script:
    """
    Rscript $projectDir/bin/tidy_bvbrc.R \
    --genome ${params.dbdir}/bvbrc/${bvbrc_db_ver}/BVBRC_genome.csv \
    --amr ${params.dbdir}/bvbrc/${bvbrc_db_ver}/BVBRC_genome_amr.csv \
    """
}

workflow {
    // prepare a tidy database of measured resistance data
    tidy_bvbrc()

    // run QC checks
    run_qc_checks()

    // calibrate carbapenem prediction
    calibrate_carbapenem_prediction(run_qc_checks.out.aci_with_qc, tidy_bvbrc.out)
}
