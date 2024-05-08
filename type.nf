#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Pipeline version
pipeline_version = "v0.13.5.9001"

// Container versions
blast_container = "ncbi/blast:2.13.0"
busco_container = "stitam/busco:5.4.4"
diamond_container = "stitam/diamond"
kaptive_container = "stitam/kaptive:2.0.3"
kraken2_container = "stitam/kraken2:dev"
mlst_container = "staphb/mlst:2.23.0"
prodigal_container = "nanozoo/prodigal:2.6.3--2769024"
prokka_container = "stitam/prokka:1.14.5"
r_container = "stitam/r-aci:0.13"

// Database versions
resfinder_db_ver = "2.0.0"
busco_db_name = "pseudomonadales"
busco_db_ver = "odb10"
kraken2_db_name = "bacteria"
kraken2_db_ver = "v230502"

// Input parameters
params.file = "assemblies.tsv"

// Output parameters
params.resdir = "results"

// Process parameters
blast_evalue = 0.00001
blast_coverage = 0.8

// Calculated process parameters
blast_coverage_perc = 100 * blast_coverage

process annotate_genome {
    container "$prokka_container"
    storeDir "$launchDir/$params.resdir/annotate_genome"

    input:
    tuple val(sample_id), path(assembly)

    output:
    path sample_id

    script:
    """
    prokka --cpus ${task.cpus} --outdir $sample_id --prefix $sample_id $assembly
    """
}

process blast_extra_genes {
    container "$blast_container"
    storeDir "$launchDir/$params.resdir/blast_extra_genes"

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_extra_genes.tsv")

    script:
    """
    blastn \
    -query $assembly \
    -subject $projectDir/data/extra_genes.fna \
    -evalue $blast_evalue \
    -out "${sample_id}_extra_genes.tsv" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore qseq sseq"
    """
}

process blast_iseqs {
    container "$blast_container"
    storeDir "$launchDir/$params.resdir/blast_iseqs"

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_iseqs.tsv")

    script:
    """
    blastn \
    -query $assembly \
    -subject $projectDir/data/insertion_sequences_LC136852.fna \
    -evalue $blast_evalue \
    -out "${sample_id}_iseqs.tsv" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore qseq sseq"
    """
}

process blast_resfinder {
    tag "$sample_id"
    container "$diamond_container"
    storeDir "$launchDir/$params.resdir/blast_resfinder"

    input:
    tuple val(sample_id), path(orfs_fna), path(orfs_faa), path(orfs_gff)

    output:
    tuple val(sample_id), path(orfs_gff), path("${sample_id}_resfinder.tsv")

    script:
    """
    diamond blastx \
    --query $orfs_fna \
    --db "${params.dbdir}/resfinder/$resfinder_db_ver/resdb.dmnd" \
    --out "${sample_id}_resfinder.tsv" \
    --sensitive \
    --block-size 8 \
    -c1 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore \
    --max-target-seqs 1 \
    --evalue $blast_evalue \
    --query-cover $blast_coverage_perc \
    --subject-cover $blast_coverage_perc \
    --min-orf 1 \
    --id 80 \
    --header \
    --threads ${task.cpus}
    """
}

process calc_assembly_stats {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/calc_assembly_stats"

    input:
    tuple val(sample_id), path(assembly)

    output:
    path "${sample_id}.tsv"

    script:
    """
    Rscript $projectDir/bin/calc_assembly_stats.R $sample_id $assembly
    """
}

process cbind_analysis_results {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/cbind_analysis_results"

    input:
    path(assembly_stats)
    path(busco)
    path(kraken2)
    path(mlst)
    path(kaptive)
    path(resgenes)
    path(crab)

    output:
    path "analysis_results.tsv"

    script:
    """
    Rscript $projectDir/bin/cbind_analysis_results.R \
      $assembly_stats \
      $busco \
      $kraken2 \
      $mlst \
      $kaptive \
      $resgenes \
      $crab \
      $launchDir/$params.file \
      $params.Rdir
    """
}

process create_empty_resfile {
    tag "$sample_id"
    storeDir "$launchDir/$params.resdir/create_empty_resfile"

    input:
    tuple val(sample_id), path(orfs)

    output:
    tuple val(sample_id), path("${sample_id}.tsv")

    script:
    """
    touch "${sample_id}.tsv"
    """
}

process filter_output {
    container "$r_container"
    storeDir "$launchDir/results"

    input:
    path results

    output:
    path "results_short.rds"
    path "results_short.tsv"
    path "results_long.rds"
    path "results_long.tsv"

    script:
    """
    Rscript $projectDir/bin/filter_output.R $pipeline_version $results $params.Rdir
    """
}

process format_extra_genes {
    tag "$sample_id"
    container "$r_container"
    storeDir "$launchDir/$params.resdir/format_extra_genes"

    input:
    tuple val(sample_id), path(extra_genes)

    output:
    tuple val(sample_id), path("${sample_id}_format_extra_genes.tsv")

    script:
    """
    Rscript $projectDir/bin/format_extra_genes.R $projectDir/R $sample_id $extra_genes $blast_coverage
    """
}

process format_iseqs {
    tag "$sample_id"
    container "$r_container"
    storeDir "$launchDir/$params.resdir/format_iseqs"

    input:
    tuple val(sample_id), path(iseqs)

    output:
    tuple val(sample_id), path("${sample_id}_format_iseqs.tsv")

    script:
    """
    Rscript $projectDir/bin/format_iseqs.R $projectDir/R $sample_id $iseqs $blast_coverage
    """
}

process format_kaptive {
    tag "$sample_id"
    container "$r_container"
    storeDir "$launchDir/$params.resdir/format_kaptive"

    input:
    tuple val(sample_id), path(k), path(o)

    output:
    path("${sample_id}_formatted.tsv")

    script:
    """
    Rscript $projectDir/bin/format_kaptive.R $sample_id $k $o $params.Rdir
    """
}

process format_mlst {
    tag "$sample_id"
    container "$r_container"
    storeDir "$launchDir/$params.resdir/format_mlst"

    input:
    tuple val(sample_id), path(mlst)

    output:
    path("${sample_id}_formatted.tsv")

    script:
    """
    Rscript $projectDir/bin/format_mlst.R $sample_id $mlst
    """
}

process format_resgenes {
    tag "$sample_id"
    container "$r_container"
    storeDir "$launchDir/$params.resdir/format_resgenes"

    input:
    tuple val(sample_id), path(orf_gff), path(resfinder)

    output:
    tuple val(sample_id),
    path("${sample_id}_format_resgenes.tsv"),
    path("${sample_id}_classes.tsv"),
    path("${sample_id}_abs.tsv"),
    path("${sample_id}_genes.tsv")

    script:
    """
    Rscript $projectDir/bin/format_resgenes.R $projectDir/R $params.dbdir $sample_id $resfinder $blast_coverage $params.Rdir $orf_gff
    """
}

process predict_carbapenem {
    tag "$sample_id"
    container "$r_container"
    storeDir "$launchDir/$params.resdir/predict_carbapenem"
    
    input:
    tuple val(sample_id), path(resgenes), path(extra_genes), path(iseqs)

    output:
    tuple val(sample_id), path("${sample_id}_predict_carbapenem.tsv"), path("${sample_id}_iseq_debug.tsv")

    script:
    """
    Rscript $projectDir/bin/predict_carbapenem.R \
    --sample_id $sample_id \
    --resgenes $resgenes \
    --extra_genes $extra_genes \
    --insertion_sequences $iseqs \
    --resgene_groups $projectDir/data/resgene_groups.csv
    """
}

process predict_fq {
    tag "$sample_id"
    container "$r_container"
    storeDir "$launchDir/$params.resdir/predict_fq"
    // errorStrategy "ignore"

    input:
    tuple val(sample_id), path(assembly), path(extra_genes)

    output:
    tuple val(sample_id), path("${sample_id}_predict_fq.tsv")

    script:
    """
    Rscript $projectDir/bin/predict_fq.R \
    $projectDir/R \
    $sample_id \
    $assembly \
    $extra_genes \
    $projectDir/data/extra_genes.faa
    """
}

process predict_k_serotype {
    tag "$sample_id"
    container "$kaptive_container"
    storeDir "$launchDir/$params.resdir/predict_k_serotype"

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_k_table.txt")

    script:
    """
    kaptive.py \
    -a $assembly \
    -k /kaptive/reference_database/Acinetobacter_baumannii_k_locus_primary_reference.gbk \
    --no_seq_out \
    --no_json \
    -o "${sample_id}_k" \
    -t ${task.cpus}
    """
}

process predict_mlst {
    tag "$sample_id"
    container "$mlst_container"
    storeDir "$launchDir/$params.resdir/predict_mlst"

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}.tsv")

    script:
    """
    mlst \
    --label $sample_id \
    --scheme abaumannii_2 $assembly \
    --threads ${task.cpus} > "${sample_id}.tsv"
    """
}

process predict_orfs {
    tag "$sample_id"
    container "$prodigal_container"
    storeDir "$launchDir/$params.resdir/predict_orf"

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id),
    path("${sample_id}_orfs.fna"),
    path("${sample_id}_orfs.faa"),
    path("${sample_id}_orfs.gff")

    script:
    """
    prodigal -q -c -m -i $assembly \
    -d "${sample_id}_orfs.fna" \
    -a "${sample_id}_orfs.faa" \
    -f gff -o "${sample_id}_orfs.gff" &
    """
}

process predict_o_serotype {
    tag "$sample_id"
    container "$kaptive_container"
    storeDir "$launchDir/$params.resdir/predict_o_serotype"

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${sample_id}_o_table.txt")

    script:
    """
    kaptive.py \
    -a $assembly \
    -k /kaptive/reference_database/Acinetobacter_baumannii_OC_locus_primary_reference.gbk \
    --no_seq_out \
    --no_json \
    -o "${sample_id}_o" \
    -t ${task.cpus}
    """
}

process predict_crab {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/predict_crab"

    input:
    tuple val(sample_id), path(resgenes), path(carbapenem), path(fq_res_mutations)

    output:
    tuple val(sample_id), path("${sample_id}_crab.tsv")

    script:
    """
    Rscript $projectDir/bin/predict_crab.R $sample_id $resgenes $carbapenem $fq_res_mutations
    """
}

process run_busco {
    tag "sample_id"
    container "$busco_container"
    containerOptions "--no-home"
    storeDir "$launchDir/$params.resdir/run_busco"

    input:
    tuple val(sample_id), path(assembly)

    output:
    path sample_id

    script:
    """
    busco \
    -i $assembly \
    -o $sample_id \
    -m geno \
    -l $params.dbdir/busco/lineages/${busco_db_name}_${busco_db_ver} \
    --cpu ${task.cpus} \
    --offline
    """
}

process run_busco_py {
    tag "$sample_id"
    container "$busco_container"
    containerOptions "--no-home"
    storeDir "$launchDir/$params.resdir/run_busco_py"

    input:
    tuple val(sample_id), path(assembly)

    output:
    path sample_id
    path "${sample_id}_busco.tsv"

    script:
    """
    python3 $projectDir/bin/run_busco.py3 \
    -n 8 \
    -g $assembly \
    -o "${sample_id}_busco.tsv" \
    -db_name $busco_db_name \
    -db_path $params.dbdir/busco \
    -db_ver $busco_db_ver
    """
}

process run_kraken2 {
    container "$kraken2_container"
    containerOptions "--no-home"
    storeDir "$launchDir/$params.resdir/run_kraken2"

    input:
    tuple val(sample_id), path(assembly)

    output:
    path("${sample_id}_kraken2.tsv")
    path("${sample_id}_kraken2_report.tsv")

    script:
    """
    python3 $projectDir/bin/run_kraken.py3 \
    -n ${task.cpus} \
    -g $assembly \
    -o "${sample_id}_kraken2.tsv" \
    -r "${sample_id}_kraken2_report.tsv" \
    -db_name $kraken2_db_name \
    -db_path $params.dbdir/kraken2 \
    -db_ver $kraken2_db_ver
    """
}

process validate_input {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/validate_input"

    output:
    path "validated_input.csv"
    path "missing_assemblies.tsv"
    path "smallfiles.tsv"
    path "bigfiles.tsv"
    
    script:
    """
    Rscript $projectDir/bin/validate_input.R $launchDir $launchDir/$params.file
    """
}

process validate_output {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/validate_output"

    input:
    path results

    output:
    path "validated_results.tsv"

    script:
    """
    Rscript $projectDir/bin/validate_output.R $launchDir/$params.file $results $pipeline_version
    """
}

workflow {
    // validate input
    validate_input()
    assembly_tuple_ch = validate_input.out[0] | splitCsv(header: false, skip: 1)
    // // annotate genomes
    // assembly_tuple_ch | annotate_genome
    // calculate assembly stats
    assembly_tuple_ch | calc_assembly_stats
    // bind assembly stats for all assemblies
    assembly_stats = calc_assembly_stats.out.collectFile(
        name: "assembly_stats.tsv",
        newLine: false,
        keepHeader: true,
        storeDir: "$launchDir/results"
    )
    // run busco
    // assembly_tuple_ch | run_busco
    assembly_tuple_ch | run_busco_py
    // bind busco data for all assemblies
    busco = run_busco_py.out[1].collectFile(
        name: "busco.tsv",
        newLine: false,
        keepHeader: true,
        storeDir: "$launchDir/results"
    )
    // run kraken2
    assembly_tuple_ch | run_kraken2
    // bind kraken2 data for all assemblies
    kraken2 = run_kraken2.out[0].collectFile(
        name: "kraken2.tsv",
        newLine: false,
        keepHeader: true,
        storeDir: "$launchDir/results"
    )
    // predict mlst
    assembly_tuple_ch | predict_mlst | format_mlst
    // bind mlst data for all assemblies
    mlst = format_mlst.out.collectFile(
        name: "mlst.tsv",
        newLine: false,
        keepHeader: true,
        storeDir: "$launchDir/results"
    )
    // predict k and o serotype
    assembly_tuple_ch | predict_k_serotype
    assembly_tuple_ch | predict_o_serotype
    predict_k_serotype.out.join(predict_o_serotype.out) | format_kaptive
    // bind kaptive data for all assemblies
    kaptive = format_kaptive.out.collectFile(
        name: "kaptive.tsv",
        newLine: false,
        keepHeader: true,
        storeDir: "$launchDir/results"
    )
    // predict orfs and resgenes
    assembly_tuple_ch | predict_orfs
    predict_orfs.out | filter{ it[1].size() > 0 } | blast_resfinder
    predict_orfs.out | filter{ it[1].size() == 0} | create_empty_resfile
    blast_resfinder.out.concat(create_empty_resfile.out) | format_resgenes
    // bind resfinder data for all assemblies
    resgenes = format_resgenes.out.map{it[1]}.collectFile(
        name: "resgenes.tsv",
        newLine: false,
        keepHeader: true,
        storeDir: "$launchDir/results"
    )
    resgenes_by_class = format_resgenes.out.map{it[2]}.collectFile(
        name: "resgenes_by_class.tsv",
        newLine: false,
        keepHeader: true,
        storeDir: "$launchDir/results"
    )
    resgenes_by_ab = format_resgenes.out.map{it[3]}.collectFile(
        name: "resgenes_by_ab.tsv",
        newLine: false,
        keepHeader: true,
        storeDir: "$launchDir/results"
    )
    resgenes_by_gene = format_resgenes.out.map{it[4]}.collectFile(
        name: "resgenes_by_gene.tsv",
        newLine: false,
        keepHeader: true,
        storeDir: "$launchDir/results"
    )
    // find extra genes and ISAba1 insertion sequences
    assembly_tuple_ch | blast_iseqs | format_iseqs
    assembly_tuple_ch | blast_extra_genes | format_extra_genes
    // predict point mutations indicative of fluoroquinolone resistance
    assembly_tuple_ch.join(format_extra_genes.out) | predict_fq
    // predict genomic features indicative of carbapenem resistance
    format_resgenes.out.map{it[0,1]}.join(format_extra_genes.out.join(format_iseqs.out)) | predict_carbapenem
    predict_carbapenem.out.map{it[2]}.collectFile(
        name: "iseq_debug.tsv",
        newLine: false,
        keepHeader: true,
        storeDir: "$launchDir/results"
    )
    // crab classification
    format_resgenes.out.map{it[0,2]}.join(predict_carbapenem.out.map{it[0,1]}.join(predict_fq.out)) | predict_crab
    // bind crab data for all assemblies
    crab = predict_crab.out.map{it[1]}.collectFile(
        name: "crab.tsv",
        newLine: false,
        keepHeader: true,
        storeDir: "$launchDir/results"
    )
    // merge prediction results into a single file and validate output
    cbind_analysis_results(
        assembly_stats,
        busco,
        kraken2,
        mlst,
        kaptive,
        resgenes_by_class,
        crab) | validate_output | filter_output
   }
