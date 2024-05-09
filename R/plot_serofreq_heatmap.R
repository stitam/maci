plot_serofreq_heatmap <-function(df, padding = 5) {
  if ("collection_year" %in% names(df) == FALSE) {
    stop("Input data frame must have a variable called 'collection_year'.")
  }
  if (any(is.na(df$collecion_year))) {
    stop("'collection_year' must be avaiable for each sample.")
  }
  
  minyear = min(df$collection_year)
  maxyear = max(df$collection_year)
  regions_serotypes_geodate <- serotype_freqs(df, group_by = "region23")
  
  index_remaining <- which(regions_serotypes_geodate$ratio_cumsum > 0.9)
  remaining_count <- length(
    unique(regions_serotypes_geodate$serotype[index_remaining])
  )
  
  regions_serotypes_geodate_hm <- regions_serotypes_geodate  %>%  
    mutate(
      regions_serotypes_geodate,
      serotype2=if_else(
        ratio_cumsum>0.9, 
        paste0("Other (", remaining_count, ")"), serotype)) %>% 
    group_by(region23, serotype2) %>% 
    summarise(percent=sum(ratio*100))
  
  regions_serotypes_geodate_wide <- regions_serotypes_geodate_hm %>% 
    pivot_wider(names_from=serotype2, values_from = percent)
  regions_serotypes_geodate_wide <- subset(
    regions_serotypes_geodate_wide,
    select = c(1, 3:ncol(regions_serotypes_geodate_wide), 2)
  )
  m <- as.matrix(regions_serotypes_geodate_wide[,-1])
  rownames(m) <- regions_serotypes_geodate_wide$region23
  cn = colnames(regions_serotypes_geodate_wide[-1])
  
  regions_serotypes_geodate_bray <- regions_serotypes_geodate  %>%  
    mutate(regions_serotypes_geodate, serotype2=serotype) %>% 
    group_by(region23, serotype2) %>% 
    summarise(count=sum(count))
  
  regions_serotypes_geodate_bray_wide <- pivot_wider(
    regions_serotypes_geodate_bray, names_from=serotype2, values_from=count)
  diversity_matrix <- as.matrix(regions_serotypes_geodate_bray_wide[,-1])
  rownames(diversity_matrix) <- regions_serotypes_geodate_bray_wide$region23
  diversity_matrix[is.na(diversity_matrix)] <- 0
  beta_diversity <- vegan::vegdist(diversity_matrix, "bray")
  
  # dendogram
  my_clustering <- hclust(dist(beta_diversity))
  
  global_prevalent <- global_prevalent(df)
  
  rest <- data.frame(serotype=colnames(regions_serotypes_geodate_wide[-1])) %>% mutate(color="#000000")
  rest <- anti_join(rest,global_prevalent,  by="serotype")
  sero_color <- rbind(global_prevalent, rest)
  sero_color <- sero_color[match(colnames(regions_serotypes_geodate_wide[-1]), sero_color$serotype),]
  
  colors = circlize::colorRamp2(
    c(1.999, 2, 4.999,5,40),
    c("gray","gray50","gray50","#DBF3FF","royalblue4")
  )
  
  # breaks for plot legend
  gbreaks <- c(0,2,5)
  maxm <- 10 * ceiling(max(m, na.rm = TRUE) / 10)
  if (maxm >= 0) {
    gbreaks <- c(gbreaks, seq(from = 10, to = maxm, by = 10))
  }
  
  index_prevalent <- which(sero_color$serotype %in% global_prevalent$serotype)
  index_prevalent <- c(index_prevalent, nrow(sero_color))
  
  xtext_linecolors <- rep(
    c("#1F1E53", "#925637"),
    times = ceiling(length(index_prevalent)/2)
  )
  xtext_linecolors <- xtext_linecolors[1:length(index_prevalent)]
  
  heatmap_region <- ComplexHeatmap::Heatmap(
    m,
    cluster_rows = my_clustering,
    cluster_columns = F,
    row_dend_width = unit(0.4, "cm"),
    row_dend_gp = grid::gpar(lwd = 0.1),
    row_names_gp = grid::gpar(fontsize = 5),
    col=colors,
    heatmap_legend_param = list(
      title="Percent",
      title_gp = grid::gpar(fontsize = 5),
      at = gbreaks,
      labels_gp = grid::gpar(fontsize = 5),
      legend_height = unit(1, "cm"),
      legend_width = unit(2, "mm"),
      color_bar = "continuous"
    ),
    rect_gp = grid::gpar(col = "grey80", lwd = 0.5),
    na_col="white",
    show_column_names = FALSE,
    row_names_side = "left",
    bottom_annotation = ComplexHeatmap::HeatmapAnnotation(
      text = ComplexHeatmap::anno_mark(
        at = index_prevalent,
        labels = sero_color$serotype[index_prevalent],
        which = "column",
        side = "bottom",
        lines_gp = grid::gpar(lwd = 0.5, col = xtext_linecolors ),
        labels_gp = grid::gpar(fontsize = 5),
        #labels_gp = grid::gpar(fontsize = 5, col=sero_color$color[index_prevalent]),
        labels_rot = 60,
        padding = unit(3, "mm")
      )
    )
  ) 
  
  g1 <- ggplotify::as.ggplot(heatmap_region) + 
    geom_label(
      aes(x = 0.82, y = 0.89, label = paste0(minyear, "-", maxyear)),
      fill = "#FFFFFF",
      size = 5/.pt,
      label.size = 0.1
    )
  
  g1$labels$label <- paste0(minyear, "-", maxyear)
  g1$layers[[3]]$constructor[[2]]$label <- paste0(minyear, "-", maxyear)
  g1$layers[[3]]$mapping$label <- paste0(minyear, "-", maxyear)
  
  return(g1)
  
}