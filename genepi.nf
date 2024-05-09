#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// TODO info about the workflow here, after typing, what output it will generate

// Container versions
r_container = "stitam/r-aci:0.14"

downsampling_strategy = Channel.of(params.downsampling_strategy)

// Input parameters
params.minyear = 2009
params.minyear_recent = 2016
params.mincount = 50

// Input files
params.assembly_summary = "${launchDir}/aci_study.rds"
params.global_rr = "input/rr_global_no_focus/relative_risks_type2.rds"
params.global_tree_nwk = "${launchDir}/data/global_ST2_tree/dated_tree.nwk"
params.global_tree_rds = "${launchDir}/data/global_ST2_tree/global_tree.rds"
params.phage_sensitivity = "${launchDir}/Ab_all_strains_phages_spotassay_PFU.tsv"
params.poppunk_clusters = "null"
params.regional_rr = "input/rr_regional/relative_risks_type3.rds"
params.trees = "${launchDir}/all_dated_trees.rds"

// Output parameters
params.resdir = "results"

// Processes 

process analyse_diversity {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/analyse_diversity"

    input:
    tuple val(strategy), path(aci_tbl)

    output:
    path "diversity_stats_ds_${strategy}.tsv"
    path "log_${strategy}.txt"

    script:
    """
    Rscript $projectDir/bin/analyse_diversity.R \
    --file $aci_tbl \
    --project_dir $projectDir \
    --downsampling_strategy $strategy \
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
      path("TableS2A_top_serotypes_region23_ds_${strategy}.tsv"), \
      path("global_or_prevalent_serotypes_region23_ds_${strategy}.tsv"), \
      path("top_serotypes_country_ds_${strategy}.tsv"), \
      path("country_comparisons_ds_${strategy}.tsv"), \
      path("serotypes_country_year_ds_${strategy}.tsv"), \
      path("meta.rds"), \
      emit: top_serotypes
    path "top_serotypes_summary_region23_ds_${strategy}.tsv"
    path "serotypes_region23_year_ds_${strategy}.tsv"
    path "global_or_prevalent_overlap_region23_ds_${strategy}.tsv"
    path "TableS2_prevalent_overlap_region23_ds_${strategy}.tsv"
    path "log.txt"

    script:
    """
    Rscript ${projectDir}/bin/calc_serotype_freqs.R \
    --project_dir $projectDir \
    --file $aci_tbl \
    --shortlist $params.geographic_locations \
    --downsampling_strategy $strategy \
    --minyear $params.minyear \
    --minyear_recent $params.minyear_recent
    """
}

process downsample_isolates {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/downsample_isolates"

    input:
    path aci_filtered
    each strategy

    output:
    tuple val(strategy), path("aci_crab_ds_${strategy}.tsv"), emit: tbl
    path "country_counts_${strategy}.tsv"

    script:
    """
    Rscript $projectDir/bin/downsample_isolates.R \
    --project_dir $projectDir \
    --aci_path $aci_filtered \
    --pp_path $params.poppunk_clusters \
    --downsampling_strategy $strategy \
    --geographic_locations $params.geographic_locations \
    --population_sampling_rate $params.population_sampling_rate
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
    --shortlist $params.geographic_locations
    """
}

process filter_crab {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/filter_crab"

    input:
    tuple val(strategy), path(aci_tbl)

    output:
    tuple val(strategy), path("aci_ds_${strategy}_crab.rds")

    script:
    """
    Rscript $projectDir/bin/filter_crab.R \
    $projectDir \
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
    storeDir "$launchDir/$params.resdir/plot_boxplot_global_serotypes/${strategy}"

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
    path "boxplot_global_serotypes_ds_${strategy}.pdf"
    path "boxplot_global_serotypes_ds_${strategy}.png"
    path "boxplot_global_serotypes_ds_${strategy}.rds", emit: boxplot_global_serotypes

    script:
    """
    Rscript $projectDir/bin/plot_boxplot_global_serotypes.R \
    --file $top_serotypes_region23 \
    --regions $params.geographic_locations \
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
    --project_dir $projectDir \
    --file $aci \
    --country_file $top_serotypes_country
    """
}

process plot_crab_over_time {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_crab_over_time"

    input:
    path aci_tbl

    output:
    path "crab_over_time.pdf"
    path "crab_over_time.png"
    path "crab_over_time.rds", emit: crab_over_time

    script:
    """
    Rscript $projectDir/bin/plot_crab_over_time.R $aci_tbl
    """
}

process plot_type_diversity {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_type_diversity/${strategy}"

    input:
    path world_map
    path crab_over_time
    tuple \
      val(strategy), \
      path(heatmap_region23), \
      path(serotop_by_region), \
      path(morisita_countries)
    
    output:
    path "Fig2.pdf"
    path "Fig2.png"

    script:
    """
    Rscript $projectDir/bin/plot_type_diversity.R \
    --fig1A $world_map \
    --fig1B $crab_over_time \
    --fig1C $serotop_by_region \
    --fig1D $heatmap_region23 \
    --fig1E $morisita_countries \
    --regions $params.geographic_locations
    """
}

process plot_temporal_dynamics {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/$params.resdir/plot_temporal_dynamics/${strategy}"

    input:
    tuple \
      val(strategy), \
      path(sero_over_time), \
      path(morisita_histograms)
    
    output:
    path "Fig3.pdf"
    path "Fig3.png"

    script:
    """
    Rscript $projectDir/bin/plot_temporal_dynamics.R \
    --figA $sero_over_time \
    --figB $morisita_histograms \
    --figC $params.global_tree_rds \
    --figD ${launchDir}/$params.global_rr \
    --figE ${launchDir}/$params.regional_rr
    """
}

process plot_phage_sensitivity {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_phage_sensitivity"

    input:
    tuple path(sensitivity_heatmap_pdf), path(sensitivity_heatmap_png)
    path sensitivity_dissimilarity
    path phylodist_phagedist

    output:
    path "Fig4.pdf"
    path "Fig4.png"

    script:
    """
    Rscript $projectDir/bin/plot_phage_sensitivity.R \
    --figA_pdf $sensitivity_heatmap_pdf \
    --figA_png $sensitivity_heatmap_png \
    --figB $sensitivity_dissimilarity \
    --figC $phylodist_phagedist
    """
}

process plot_genome_metrics {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_genome_metrics"

    input:
    path aci_tbl

    output:
    path "genome_metrics.pdf"
    path "genome_metrics.png"

    script:
    """
    Rscript $projectDir/bin/plot_genome_metrics.R --file $aci_tbl
    """
}

process plot_heatmap {
    container "$r_container"
    containerOptions "--no-home"
    storeDir "$launchDir/$params.resdir/plot_heatmap/${strategy}"

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
    path "heatmap_region23_ds_${strategy}.pdf"
    path "heatmap_region23_ds_${strategy}.png"
    tuple val(strategy), path("heatmap_region23_ds_${strategy}.rds"), emit: heatmap_region23
    path "heatmap_region23_vertical_ds_${strategy}.pdf"
    path "heatmap_region23_vertical_ds_${strategy}.png"
    path "heatmap_region23_vertical_ds_${strategy}.rds", emit: heatmap_region23_vertical
    path "heatmap_country_ds_${strategy}.pdf"
    path "heatmap_country_ds_${strategy}.png"
    path "heatmap_country_ds_${strategy}.rds", emit: heatmap_country

    script:
    """
    Rscript $projectDir/bin/plot_heatmap.R \
    --region_file $top_serotypes_region23 \
    --serotype_file $global_prevalent \
    --country_file $top_serotypes_country \
    --shortlist $params.geographic_locations \
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
    path "morisita_country_year_ds_${strategy}_crab.tsv"
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
    storeDir "$launchDir/$params.resdir/plot_overlap_morisita/${strategy}"

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
    path "Morisita_countries.pdf"
    path "Morisita_countries.png"
    path "Morisita_countries.tsv"
    tuple val(strategy), path("Morisita_countries.rds"), emit: morisita_countries

    script:
    """
    Rscript $projectDir/bin/plot_overlap_morisita.R \
    --file_path $country_comparisons
    """
}

process plot_phylodist_phagedist {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_phylodist_phagedist"

    input:
    path aci_tbl

    output:
    path "phylodist_phagedist.pdf"
    path "phylodist_phagedist.png"
    path "phylodist_phagedist.rds", emit: phylodist_phagedist
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
    storeDir "$launchDir/$params.resdir/plot_prevalence_piechart/${strategy}"

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
    path "prevalence_piechart.pdf"
    path "prevalence_piechart.png"
    tuple val(strategy), path("prevalence_piechart.rds"), emit: prevalence_piechart

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
    storeDir "$launchDir/$params.resdir/plot_rarefaction_curve/${strategy}"

    input:
    tuple val(strategy), path(aci_tbl)

    output:
    path "rarecurve_crab_ds_${strategy}_minyear_${params.minyear_recent}_regions.pdf"
    path "rarecurve_crab_ds_${strategy}_minyear_${params.minyear_recent}_regions.png"
    path "rarecurve_crab_ds_${strategy}_minyear_${params.minyear_recent}_regions.rds"
    path "rarecurve_Europe_crab_ds_${strategy}_minyear_${params.minyear_recent}_countries.pdf"
    path "rarecurve_Europe_crab_ds_${strategy}_minyear_${params.minyear_recent}_countries.png"
    path "rarecurve_not_Europe_crab_ds_${strategy}_minyear_${params.minyear_recent}_countries.pdf"
    path "rarecurve_not_Europe_crab_ds_${strategy}_minyear_${params.minyear_recent}_countries.png"
    path "log.txt"

    script:
    """
    Rscript $projectDir/bin/plot_rarefaction_curve.R \
    --project_dir $projectDir \
    --file $aci_tbl \
    --minyear $params.minyear_recent \
    --regions $params.geographic_locations
    """
}

process plot_sensitivity_dissimilarity {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_sensitivity_dissimilarity"

    output:
    path "sensitivity_dissimilarity.pdf"
    path "sensitivity_dissimilarity.png"
    path "sensitivity_dissimilarity.rds", emit: sensitivity_dissimilarity
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
    storeDir "$launchDir/$params.resdir/plot_sensitivity_heatmap"

    output:
    path "heatmap_phages_PFU_composed_country_order_noLfw_cities_old_coloring.pdf"
    path "heatmap_phages_PFU_composed_country_order_noLfw_cities_old_coloring.png"
    tuple \
      path("heatmap_phages_PFU_composed_country_order_noLfw_cities_old_coloring_pdf.rds"),
      path("heatmap_phages_PFU_composed_country_order_noLfw_cities_old_coloring_png.rds"),
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
    storeDir "$launchDir/$params.resdir/plot_sero_over_time/${strategy}"

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
    path "sero_over_time_ds_${strategy}.pdf"
    path "sero_over_time_ds_${strategy}.png"
    tuple val(strategy), path("sero_over_time_ds_${strategy}.rds"), emit: sero_over_time
    path "log_${strategy}.txt"

    script:
    """
    Rscript $projectDir/bin/plot_sero_over_time.R \
    --project_dir $projectDir \
    --file $aci_tbl \
    --serotype_file $global_prevalent \
    --minyear $params.minyear \
    --maxyear $params.maxyear
    """
}

process plot_serotop {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_serotop/${strategy}"

    input:
    tuple \
      val(strategy), \
      path(serotop_continent), \
      path(serotop_region23), \
      path(serotop_country)
    
    output:
    path "serotop_region23_ds_${strategy}.pdf"
    path "serotop_region23_ds_${strategy}.png"
    tuple val(strategy), path("serotop_region23_ds_${strategy}.rds"), emit: serotop_by_region
    path "serotop_country_ds_${strategy}.pdf"
    path "serotop_country_ds_${strategy}.png"
    path "serotop_country_ds_${strategy}.rds"

    script:
    """
    Rscript $projectDir/bin/plot_serotop.R \
    --serotop_region23 $serotop_region23 \
    --serotop_country $serotop_country \
    --regions $params.geographic_locations
    """
}

process plot_transmission_lineages {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_transmission_lineages"

    output:
    path "transmission_lineages.pdf"
    path "transmission_lineages.png"
    path "transmission_lineages.rds"
    path "chain_length_from_closest_transmission.tsv"
    path "longest_chain_from_closest_transmission_no_singletons.tsv"

    script:
    """
    Rscript $projectDir/bin/plot_transmission_lineages.R \
    --tree $params.global_tree_nwk \
    --tree_tbl $params.pastml_tbl
    """
}


process plot_world_map {
    container "$r_container"
    storeDir "$launchDir/$params.resdir/plot_world_map"  

    input:
    path aci_tbl

    output:
    path "world_map.pdf"
    path "world_map.png"
    path "world_map.rds", emit: world_map

    script:
    """
    Rscript $projectDir/bin/plot_world_map.R \
    --file $aci_tbl \
    --regions $params.geographic_locations
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
      path("serotop_continent_ds_${strategy}.tsv"), \
      path("serotop_region23_ds_${strategy}.tsv"), \
      path("serotop_country_ds_${strategy}.tsv"), \
      emit: serotop_tbls
    path "log_${strategy}.txt"

    script:
    """
    Rscript ${projectDir}/bin/prep_serotop_tbls.R \
    --project_dir $projectDir \
    --file $aci_tbl \
    --strategy $strategy \
    --minyear $params.minyear_recent
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
    --file $params.assembly_summary \
    --taxid $params.taxid
    """
}


workflow {
    // run QC checks
    run_qc_checks()

    // FILTER, DOWNSAMPLE, PREPARE TABLES

    // add new variable to flag filtered isolates (by assembly QC and metadata)
    filter_assemblies(run_qc_checks.out.aci_with_qc)

    // add new variable to flag downsampled crab isolates
    downsample_isolates(filter_assemblies.out.aci_filtered, downsampling_strategy)
    // // analyse diversity
    // downsample_isolates.out.tbl | analyse_diversity

    // serotype frequencies and cumulative serotype frequencies for crab isolates
    downsample_isolates.out.tbl | calc_serotype_freqs
    // // logistic regressions for crab isolates
    logistic_regressions(
        downsample_isolates.out.tbl.join(calc_serotype_freqs.out.top_serotypes)
    )
    // morisita indices for selected countries and selected timeframes
    calc_serotype_freqs.out.top_serotypes | plot_morisita_histograms

    // PLOTS FROM FILTERED DATASET

    // plot world map
    filter_assemblies.out.aci_filtered | plot_world_map
    // // prep qc plots - compare NCBI and current study
    filter_assemblies.out.aci_filtered | plot_genome_metrics
    // plot crab/non-crab counts over time
    filter_assemblies.out.aci_filtered | plot_crab_over_time

    // // PLOTS FROM FILTERED, CRAB, DOWNSAMPLED DATASET

    // plot serotype frequencies over time
    downsample_isolates.out.tbl.join(calc_serotype_freqs.out.top_serotypes) | plot_sero_over_time

    // serotop tbl and plot for crab isolates
    downsample_isolates.out.tbl | prep_serotop_tbls
    plot_serotop(prep_serotop_tbls.out.serotop_tbls)

    // heatmap for crab isolates
    plot_heatmap(calc_serotype_freqs.out.top_serotypes)

    // boxplots for global crab serotypes
    plot_boxplot_global_serotypes(calc_serotype_freqs.out.top_serotypes)

    // rarefaction curves for crab isolates
    downsample_isolates.out.tbl | plot_rarefaction_curve

    // // boxplots of bray-curtis distances
    // input_bray_curtis = downsample_isolates.out.tbl.join(calc_serotype_freqs.out.top_serotypes)

    // input_bray_curtis | plot_bray_curtis

    // boxplots of 1. mean prevalence of overlapping serotypes 2. morisita indices
    plot_overlap_morisita(calc_serotype_freqs.out.top_serotypes)

    // Fig1

    fig1_sub = plot_heatmap.out.heatmap_region23.join(
        plot_serotop.out.serotop_by_region.join(
            plot_overlap_morisita.out.morisita_countries
        )
    )

    plot_type_diversity(
        plot_world_map.out.world_map,
        plot_crab_over_time.out.crab_over_time,
        fig1_sub
    )

    plot_prevalence_piechart(
        calc_serotype_freqs.out.top_serotypes.join(
            plot_sero_over_time.out.sero_over_time
        )
    )

    plot_transmission_lineages()

    plot_temporal_dynamics(
        plot_sero_over_time.out.sero_over_time.join(
            plot_morisita_histograms.out.morisita_histograms
        )
    )

    plot_sensitivity_heatmap()

    filter_assemblies.out.aci_filtered | plot_phylodist_phagedist

    plot_sensitivity_dissimilarity() 

    plot_phage_sensitivity(
        plot_sensitivity_heatmap.out.sensitivity_heatmap,
        plot_sensitivity_dissimilarity.out.sensitivity_dissimilarity,
        plot_phylodist_phagedist.out.phylodist_phagedist
    )

}