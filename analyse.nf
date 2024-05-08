#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// This pipeline includes processes that are performed after a unified aci
// prediction results table has been created.

// Container versions
r_container = "stitam/r-aci:0.13"
python_container = "mesti90/hgttree:2.11"

collapse_strategy = Channel.of("none", "geodate", "poppunk")

// Input parameters
params.assembly_summary = "${launchDir}/aci_study.rds"
params.poppunk_clusters = "${launchDir}/ppdb_clusters.csv"
params.phage_sensitivity = "${launchDir}/Ab_all_strains_phages_spotassay_PFU.tsv"
params.global_tree = "${launchDir}/input/global_ST2_tree/dated_tree.rds"
params.trees = "${launchDir}/all_dated_trees.rds"
params.minyear = 2009
params.minyear_recent = 2016
params.mincount = 50

// Input files
params.global_rr = "input/rr_global_no_focus/relative_risks_type2.rds"
params.regional_rr = "input/rr_regional/relative_risks_type3.rds"

// Output parameters
params.resdir = "results"

// Removed from workflow
process add_measured_resistance {
    container "$python_container"
    containerOptions "--no-home"
    storeDir "$launchDir/$params.resdir/add_measured_resistance"

    input:
    path aci_filtered

    output:
    path "aci_with_bvbrc.tsv", emit: aci_with_bvbrc

    script:
    """
    python3 $projectDir/bin/add_measured_resistance.py3 \
    -bvbrc_genome ${params.dbdir}/bvbrc/${bvbrc_db_ver}/BVBRC_genome.csv \
    -bvbrc_amr ${params.dbdir}/bvbrc/${bvbrc_db_ver}/BVBRC_genome_amr.csv \
    -assemblies $aci_filtered \
    -output aci_with_bvbrc.tsv
    """
}

process analyse_diversity {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/analyse_diversity"

    input:
    tuple val(strategy), path(aci_tbl)

    output:
    path "diversity_stats_collapse_${strategy}.tsv"
    path "log_${strategy}.txt"

    script:
    """
    Rscript $projectDir/bin/analyse_diversity.R \
    --file $aci_tbl \
    --project_dir $projectDir \
    --collapse_strategy $strategy \
    --minyear $params.minyear \
    --minyear_recent $params.minyear_recent
    """
}

process calc_serotype_freqs {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/calc_serotype_freqs/${strategy}"

    input:
    tuple val(strategy), path(aci_tbl)

    output:
    tuple \
      val(strategy), \
      path("TableS2A_top_serotypes_region23_collapse_${strategy}.tsv"), \
      path("global_or_prevalent_serotypes_region23_collapse_${strategy}.tsv"), \
      path("top_serotypes_country_collapse_${strategy}.tsv"), \
      path("country_comparisons_collapse_${strategy}.tsv"), \
      path("serotypes_country_year_collapse_${strategy}.tsv"), \
      path("meta.rds"), \
      emit: top_serotypes
    path "top_serotypes_summary_region23_collapse_${strategy}.tsv"
    path "serotypes_region23_year_collapse_${strategy}.tsv"
    path "global_or_prevalent_overlap_region23_collapse_${strategy}.tsv"
    path "TableS2_prevalent_overlap_region23_collapse_${strategy}.tsv"
    path "log.txt"

    script:
    """
    Rscript ${projectDir}/bin/calc_serotype_freqs.R \
    --file $aci_tbl \
    --shortlist ${projectDir}/data/geographic_locations_in_study.tsv \
    --collapse_strategy $strategy \
    --minyear $params.minyear \
    --minyear_recent $params.minyear_recent
    """
}

process collapse_outbreaks {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/collapse_outbreaks"

    input:
    path aci_filtered
    each strategy

    output:
    tuple val(strategy), path("aci_collapse_${strategy}.rds")

    script:
    """
    Rscript $projectDir/bin/collapse_outbreaks.R \
    --project_dir $projectDir \
    --aci_path $aci_filtered \
    --pp_path $params.poppunk_clusters \
    --collapse_strategy $strategy
    """
}

process describe_poppunk {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/describe_poppunk"

    input:
    path aci_filtered

    output:
    path "aci_summary.rds"
    path "stats_by_pp.rds"
    path "stats_by_pp.tsv"
    path "poppunk_cluster_sizes.png"
    path "number_of_K_serotypes_within_a_poppunk_cluster.png"
    path "number_of_countries_within_a_poppunk_cluster.png"


    script:
    """
    Rscript $projectDir/bin/describe_poppunk.R \
    $aci_filtered \
    $params.poppunk_clusters
    """
}

process filter_assemblies {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/filter_assemblies"

    input:
    path aci_with_qc

    output:
    path "aci_filtered.rds", emit: aci_filtered
    path "aci_filtered.tsv", emit: aci_filtered_tsv
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/filter_assemblies.R \
    --project_dir $projectDir \
    --file $aci_with_qc \
    --shortlist $projectDir/data/geographic_locations_in_study.tsv
    """
}

process filter_crab {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/filter_crab"

    input:
    tuple val(strategy), path(aci_tbl)

    output:
    tuple val(strategy), path("aci_collapse_${strategy}_crab.rds")

    script:
    """
    Rscript $projectDir/bin/filter_crab.R \
    $aci_tbl \
    $strategy
    """
}

process logistic_regressions {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/logistic_regressions/$strategy"

    input:
    tuple \
      val(strategy), \
      path(aci_tbl), \
      path(C), \
      path(global_prevalent), \
      path(E), \
      path(F), \
      path(G), \
      path(H)

    output:
    path "figures"
    path "logistic_regression_results.rds"
    path "logistic_regression_results_tidy.rds"
    path "logistic_regression_results_tidy.tsv"
    path "logistic_regression_results.tsv"
    path "logistic_regression_results_with_fdr.rds"
    path "TableS3_logistic_regression_results_with_fdr.tsv"

    script:
    """
    Rscript $projectDir/bin/logistic_regressions.R \
    --project_dir $projectDir \
    --aci_path $aci_tbl \
    --minyear_threshold $params.minyear \
    --mincount_threshold $params.mincount \
    --selected_serotypes $global_prevalent
    """
}

process plot_boxplot_global_serotypes {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_FigS3_boxplot_global_serotypes/${strategy}"

    input:
    tuple \
      val(strategy), \
      path(top_serotypes_region23), \
      path(global_prevalent), \
      path(top_serotypes_country), \
      path(E), \
      path(F), \
      path(G)

    output:
    path "FigS3_boxplot_global_serotypes_collapse_${strategy}.pdf"
    path "FigS3_boxplot_global_serotypes_collapse_${strategy}.png"
    path "FigS3_boxplot_global_serotypes_collapse_${strategy}.rds", emit: boxplot_global_serotypes

    script:
    """
    Rscript $projectDir/bin/plot_boxplot_global_serotypes.R \
    --file $top_serotypes_region23 \
    --regions $projectDir/data/geographic_locations_in_study.tsv \
    --serotype_file $global_prevalent
    """
}

process plot_bray_curtis {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_bray_curtis/$strategy"

    input:
    tuple \
      val(strategy), \
      path(aci), \
      path(C), \
      path(D), \
      path(top_serotypes_country), \
      path(F), \
      path(G), \
      path(H)

    output:
    path "bray_curtis.pdf"
    path "bray_curtis.png"
    path "bray_curtis.rds", emit: bray_curtis

    script:
    """
    Rscript $projectDir/bin/plot_bray_curtis.R \
    --file $aci \
    --country_file $top_serotypes_country
    """
}

process plot_crab_over_time {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_Fig1B_crab_over_time"

    input:
    path aci_tbl

    output:
    path "Fig1B_crab_over_time.pdf"
    path "Fig1B_crab_over_time.png"
    path "Fig1B_crab_over_time.rds", emit: crab_over_time

    script:
    """
    Rscript $projectDir/bin/plot_crab_over_time.R $aci_tbl
    """
}

process plot_fig1 {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_Fig1/${strategy}"

    input:
    path world_map
    path crab_over_time
    tuple \
      val(strategy), \
      path(heatmap_region23), \
      path(serotop_by_region), \
      path(morisita_countries)
    
    output:
    path "Fig1.pdf"
    path "Fig1.png"

    script:
    """
    Rscript $projectDir/bin/plot_fig1.R \
    --fig1A $world_map \
    --fig1B $crab_over_time \
    --fig1C $serotop_by_region \
    --fig1D $heatmap_region23 \
    --fig1E $morisita_countries \
    --regions $projectDir/data/geographic_locations_in_study.tsv
    """
}

process plot_fig2 {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/$params.resdir/plot_Fig2/${strategy}"

    input:
    tuple \
      val(strategy), \
      path(sero_over_time), \
      path(global_tree), \
      path(morisita_histograms)
    
    output:
    path "Fig2.pdf"
    path "Fig2.png"

    script:
    """
    Rscript $projectDir/bin/plot_fig2.R \
    --figA $sero_over_time \
    --figB $morisita_histograms \
    --figC $global_tree \
    --figD ${launchDir}/$params.global_rr \
    --figE ${launchDir}/$params.regional_rr
    """
}

process plot_fig3 {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_Fig3"

    input:
    tuple path(sensitivity_heatmap_pdf), path(sensitivity_heatmap_png)
    path sensitivity_dissimilarity
    path phylodist_phagedist

    output:
    path "Fig3.pdf"
    path "Fig3.png"

    script:
    """
    Rscript $projectDir/bin/plot_fig3.R \
    --figA_pdf $sensitivity_heatmap_pdf \
    --figA_png $sensitivity_heatmap_png \
    --figB $sensitivity_dissimilarity \
    --figC $phylodist_phagedist
    """
}

process plot_genome_metrics {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_FigS1_genome_metrics"

    input:
    path aci_tbl

    output:
    path "FigS1_genome_metrics.pdf"
    path "FigS1_genome_metrics.png"

    script:
    """
    Rscript $projectDir/bin/plot_genome_metrics.R --file $aci_tbl
    """
}

process plot_global_tree {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/$params.resdir/plot_global_tree/${strategy}"
    
    input:
    tuple val(strategy), path(sero_over_time)
    each assemblies_filtered

    output:
    path "global_tree.pdf"
    path "global_tree.png"
    tuple val(strategy), path("global_tree.rds"), emit: global_tree

    script:
    """
    Rscript $projectDir/bin/plot_global_tree.R \
    --project_dir $projectDir \
    --tree $params.global_tree \
    --metadata $assemblies_filtered \
    --sample_size 0    \
    --sero_over_time $sero_over_time \
    --regions $projectDir/data/geographic_locations_in_study.tsv
    """
}

process plot_heatmap {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/$params.resdir/plot_Fig1D_S2C_heatmap/${strategy}"

    input:
    tuple \
      val(strategy), \
      path(top_serotypes_region23), \
      path(global_prevalent), \
      path(top_serotypes_country), \
      path(E), \
      path(F), \
      path(metadata)

    output:
    path "Fig1D_heatmap_region23_collapse_${strategy}.pdf"
    path "Fig1D_heatmap_region23_collapse_${strategy}.png"
    tuple val(strategy), path("Fig1D_heatmap_region23_collapse_${strategy}.rds"), emit: heatmap_region23
    path "Fig1D_heatmap_region23_vertical_collapse_${strategy}.pdf"
    path "Fig1D_heatmap_region23_vertical_collapse_${strategy}.png"
    path "Fig1D_heatmap_region23_vertical_collapse_${strategy}.rds", emit: heatmap_region23_vertical
    path "FigS2C_heatmap_country_collapse_${strategy}.pdf"
    path "FigS2C_heatmap_country_collapse_${strategy}.png"
    path "FigS2C_heatmap_country_collapse_${strategy}.rds", emit: heatmap_country

    script:
    """
    Rscript $projectDir/bin/plot_heatmap.R \
    --region_file $top_serotypes_region23 \
    --serotype_file $global_prevalent \
    --country_file $top_serotypes_country \
    --shortlist $projectDir/data/geographic_locations_in_study.tsv \
    --metadata $metadata
    """
}

process plot_morisita_histograms {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_morisita_histograms"

    input:
    tuple \
      val(strategy), \
      path(top_serotypes_region23), \
      path(global_prevalent), \
      path(top_serotypes_country), \
      path(comparisons), \
      path(serotype_counts_per_year), \
      path(G)

    output:
    path "morisita_country_year_collapse_${strategy}_crab.tsv"
    path "Morisita_histograms_${strategy}_crab.pdf"
    path "Morisita_histograms_${strategy}_crab.png"
    tuple val(strategy), path("Morisita_histograms_${strategy}_crab.rds"), emit: morisita_histograms

    script:
    """
    Rscript $projectDir/bin/plot_morisita_histograms.R \
    --all_counts $top_serotypes_country \
    --by_year $serotype_counts_per_year \
    --comparisons $comparisons \
    --minyear $params.minyear \
    --minyear_recent $params.minyear_recent
    """
}

process plot_overlap_morisita {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_Fig1E_overlap_morisita/${strategy}"

    input:
    tuple \
      val(strategy), \
      path(top_serotypes_region23), \
      path(global_prevalent), \
      path(top_serotypes_country), \
      path(country_comparisons), \
      path(F), \
      path(G)

    output:
    path "Mean_prevalence_overlapping_serotypes_b.pdf"
    path "Mean_prevalence_overlapping_serotypes_b.png"
    path "Mean_prevalence_overlapping_serotypes_b.rds"
    path "Fig1E_Morisita_countries.pdf"
    path "Fig1E_Morisita_countries.png"
    path "Fig1E_Morisita_countries.tsv"
    tuple val(strategy), path("Fig1E_Morisita_countries.rds"), emit: morisita_countries

    script:
    """
    Rscript $projectDir/bin/plot_overlap_morisita.R \
    --file_path $country_comparisons
    """
}

process plot_phylodist_phagedist {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_Fig3C_phylodist_phagedist"

    input:
    path aci_tbl

    output:
    path "Fig3C_phylodist_phagedist.pdf"
    path "Fig3C_phylodist_phagedist.png"
    path "Fig3C_phylodist_phagedist.rds", emit: phylodist_phagedist
    path "phylodist_phagedist.tsv"
    path "phage_host_resistance_profiles.tsv"

    script:
    """
    Rscript $projectDir/bin/plot_phylodist_phagedist.R \
    --assemblies $aci_tbl \
    --trees $params.trees \
    --sensitivity $params.phage_sensitivity
    """
}

process plot_prevalence_piechart {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_FigS2A2_prevalence_piechart/${strategy}"

    input:
    tuple \
      val(strategy), \
      path(top_serotypes_region23), \
      path(C), \
      path(D), \
      path(E), \
      path(F), \
      path(G), \
      path(sero_over_time)

    output:
    path "FigS2A2_prevalence_piechart.pdf"
    path "FigS2A2_prevalence_piechart.png"
    tuple val(strategy), path("FigS2A2_prevalence_piechart.rds"), emit: prevalence_piechart

    script:
    """
    Rscript $projectDir/bin/plot_prevalence_piechart.R \
    --region_file $top_serotypes_region23 \
    --sero_over_time $sero_over_time \
    --minyear_recent $params.minyear_recent
    """
}

process plot_rarefaction_curve {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_FigS2A_B_rarefaction_curve/${strategy}"

    input:
    tuple val(strategy), path(aci_tbl)

    output:
    path "FigS2A_rarecurve_crab_collapse_${strategy}_minyear_${params.minyear_recent}_regions.pdf"
    path "FigS2A_rarecurve_crab_collapse_${strategy}_minyear_${params.minyear_recent}_regions.png"
    path "FigS2A_rarecurve_crab_collapse_${strategy}_minyear_${params.minyear_recent}_regions.rds"
    path "FigS2B_rarecurve_Europe_crab_collapse_${strategy}_minyear_${params.minyear_recent}_countries.pdf"
    path "FigS2B_rarecurve_Europe_crab_collapse_${strategy}_minyear_${params.minyear_recent}_countries.png"
    path "FigS2B2_rarecurve_not_Europe_crab_collapse_${strategy}_minyear_${params.minyear_recent}_countries.pdf"
    path "FigS2B2_rarecurve_not_Europe_crab_collapse_${strategy}_minyear_${params.minyear_recent}_countries.png"
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/plot_rarefaction_curve.R \
    --file $aci_tbl \
    --minyear $params.minyear_recent \
    --regions $projectDir/data/geographic_locations_in_study.tsv
    """
}

process plot_sensitivity_dissimilarity {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_Fig3B_sensitivity_dissimilarity"

    output:
    path "Fig3B_sensitivity_dissimilarity.pdf"
    path "Fig3B_sensitivity_dissimilarity.png"
    path "Fig3B_sensitivity_dissimilarity.rds", emit: sensitivity_dissimilarity
    path "histogram_3regions_bootstrap.pdf"
    
    script:
    """
    Rscript $projectDir/bin/jaccard_final.R \
    --file $params.assembly_summary \
    --sensitivity $params.phage_sensitivity
    """
}

process plot_sensitivity_heatmap {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_Fig3A_sensitivity_heatmap"

    output:
    path "Fig3A_heatmap_phages_PFU_composed_country_order_noLfw_cities_old_coloring.pdf"
    path "Fig3A_heatmap_phages_PFU_composed_country_order_noLfw_cities_old_coloring.png"
    tuple \
      path("Fig3A_heatmap_phages_PFU_composed_country_order_noLfw_cities_old_coloring_pdf.rds"),
      path("Fig3A_heatmap_phages_PFU_composed_country_order_noLfw_cities_old_coloring_png.rds"),
      emit: sensitivity_heatmap

    script:
    """
    Rscript $projectDir/bin/heatmap_PFU_5countries_final.R \
    --file $params.assembly_summary \
    --sensitivity $params.phage_sensitivity
    """
}

process plot_sero_over_time {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_Fig2A_sero_over_time/${strategy}"

    input:
    tuple \
      val(strategy), \
      path(aci_tbl), \
      path(C), \
      path(global_prevalent), \
      path(E), \
      path(F), \
      path(G), \
      path(H)

    output:
    path "Fig2A_sero_over_time_collapse_${strategy}.pdf"
    path "Fig2A_sero_over_time_collapse_${strategy}.png"
    tuple val(strategy), path("Fig2A_sero_over_time_collapse_${strategy}.rds"), emit: sero_over_time
    path "log_${strategy}.txt"

    script:
    """
    Rscript $projectDir/bin/plot_sero_over_time.R \
    --project_dir $projectDir \
    --file $aci_tbl \
    --serotype_file $global_prevalent \
    --minyear $params.minyear
    """
}

process plot_serotop {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_Fig1C_SX_serotop/${strategy}"

    input:
    tuple \
      val(strategy), \
      path(serotop_continent), \
      path(serotop_region23), \
      path(serotop_country)
    
    output:
    path "Fig1C_serotop_region23_collapse_${strategy}.pdf"
    path "Fig1C_serotop_region23_collapse_${strategy}.png"
    tuple val(strategy), path("Fig1C_serotop_region23_collapse_${strategy}.rds"), emit: serotop_by_region
    path "FigSX_serotop_country_collapse_${strategy}.pdf"
    path "FigSX_serotop_country_collapse_${strategy}.png"
    path "FigSX_serotop_country_collapse_${strategy}.rds"

    script:
    """
    Rscript $projectDir/bin/plot_serotop.R \
    --serotop_region23 $serotop_region23 \
    --serotop_country $serotop_country \
    --regions $projectDir/data/geographic_locations_in_study.tsv
    """
}

process plot_eu_trees {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_FigS6_eu_trees"

    input:
    path assemblies
    path trees
    tuple path(fig3A), path(B)
    path fig3C

    output:
    path "serotype_trees"
    path "FigS6_serotype_trees.pdf"
    path "FigS6_serotype_trees.png"
    path "tree_stats.tsv"

    script:
    """
    Rscript $projectDir/bin/plot_eu_trees.R \
    --project_dir $projectDir \
    --file $assemblies \
    --trees $trees \
    --sensitivity $params.phage_sensitivity \
    --fig3C $fig3C
    """
}


process plot_world_map {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_Fig1A_world_map"  

    input:
    path aci_tbl

    output:
    path "Fig1A_world_map.pdf"
    path "Fig1A_world_map.png"
    path "Fig1A_world_map.rds", emit: world_map

    script:
    """
    Rscript $projectDir/bin/plot_world_map.R \
    --file $aci_tbl \
    --regions $projectDir/data/geographic_locations_in_study.tsv
    """
}

process prep_serotop_tbls {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/prep_serotop_tbls"

    input:
    tuple val(strategy), path(aci_tbl)

    output:
    tuple \
      val(strategy), \
      path("serotop_continent_collapse_${strategy}.tsv"), \
      path("serotop_region23_collapse_${strategy}.tsv"), \
      path("serotop_country_collapse_${strategy}.tsv"), \
      emit: serotop_tbls
    path "log_${strategy}.txt"

    script:
    """
    Rscript ${projectDir}/bin/prep_serotop_tbls.R \
    $aci_tbl \
    $strategy \
    $params.minyear_recent
    """
}

process run_qc_checks {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/run_qc_checks"

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

workflow {
    // run QC checks
    run_qc_checks()

    // FILTER, COLLAPSE, PREPARE TABLES

    // filter data set, add new variables
    filter_assemblies(run_qc_checks.out.aci_with_qc)

    // collapse outbreaks and create new channel
    collapse_outbreaks(filter_assemblies.out.aci_filtered, collapse_strategy)
    // analyse diversity
    collapse_outbreaks.out | analyse_diversity
    // filter to crab isolates
    collapse_outbreaks.out | filter_crab
    // describe poppunk clusters
    filter_assemblies.out.aci_filtered | describe_poppunk
    // serotype frequencies and cumulative serotype frequencies for crab isolates
    filter_crab.out | calc_serotype_freqs
    // logistic regressions for crab isolates
    logistic_regressions(
        filter_crab.out.join(calc_serotype_freqs.out.top_serotypes)
    )
    // morisita indices for selected countries and selected timeframes
    calc_serotype_freqs.out.top_serotypes | plot_morisita_histograms

    // PLOTS FROM FILTERED DATASET

    // plot world map
    filter_assemblies.out.aci_filtered | plot_world_map
    // prep qc plots - compare NCBI and current study
    filter_assemblies.out.aci_filtered | plot_genome_metrics
    // plot crab/non-crab counts over time
    filter_assemblies.out.aci_filtered | plot_crab_over_time

    // PLOTS FROM FILTERED, COLLAPSED, CRAB DATASET

    // plot serotype frequencies over time
    filter_crab.out.join(calc_serotype_freqs.out.top_serotypes) | plot_sero_over_time

    // serotop tbl and plot for crab isolates
    filter_crab.out | prep_serotop_tbls
    plot_serotop(prep_serotop_tbls.out.serotop_tbls)

    // heatmap for crab isolates
    plot_heatmap(calc_serotype_freqs.out.top_serotypes)

    // boxplots for global crab serotypes
    plot_boxplot_global_serotypes(calc_serotype_freqs.out.top_serotypes)

    // rarefaction curves for crab isolates
    filter_crab.out | plot_rarefaction_curve

    // boxplots of bray-curtis distances
    input_bray_curtis = filter_crab.out.join(calc_serotype_freqs.out.top_serotypes)

    input_bray_curtis | plot_bray_curtis

    // boxplots of 1. mean prevalence of overlapping serotypes 2. morisita indices
    plot_overlap_morisita(calc_serotype_freqs.out.top_serotypes)

    // Fig1

    fig1_sub = plot_heatmap.out.heatmap_region23.join(
        plot_serotop.out.serotop_by_region.join(
            plot_overlap_morisita.out.morisita_countries
        )
    )

    plot_fig1(
        plot_world_map.out.world_map,
        plot_crab_over_time.out.crab_over_time,
        fig1_sub
    )

    plot_global_tree(
        plot_sero_over_time.out.sero_over_time,
        filter_assemblies.out.aci_filtered
    )

    plot_fig2(
        plot_sero_over_time.out.sero_over_time.join(
            plot_global_tree.out.global_tree.join(
                plot_morisita_histograms.out.morisita_histograms
            )
        )
    )
    
    plot_sensitivity_heatmap()

    plot_sensitivity_dissimilarity()

    filter_assemblies.out.aci_filtered | plot_phylodist_phagedist

    plot_fig3(
        plot_sensitivity_heatmap.out.sensitivity_heatmap,
        plot_sensitivity_dissimilarity.out.sensitivity_dissimilarity,
        plot_phylodist_phagedist.out.phylodist_phagedist
    )

    plot_prevalence_piechart(
        calc_serotype_freqs.out.top_serotypes.join(
            plot_sero_over_time.out.sero_over_time
        )
    )

    plot_eu_trees(
        filter_assemblies.out.aci_filtered,
        params.trees,
        plot_sensitivity_heatmap.out.sensitivity_heatmap,
        plot_phylodist_phagedist.out.phylodist_phagedist
    )
}