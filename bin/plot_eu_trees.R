# This script is one of multiple scripts that plot phylogenetic trees. This
# script creates a phylogenetic tree plot with fan layout. The tree will also 
# contain heatmaps for MLST, K type, country and city and will highlight tips
# that were included in phage-host laboratory tests.
#
# The script takes
# - a dated phylogenetic tree which was created by the Nextflow pipeline and
#   which includes any duplicates that were removed during tree building
# - a metadata table which contains typing results and geographical locations
#   for aci strains
# - an rds file which contains a heatmap abour phage-host laboratory test
#   results
# 
# It is expected that this script will be run after the Nextflow pipeline which
# creates the dated phylogenetic trees. The working directory for running this
# script should be the directory from which the Nextflow pipeline that build the
# phylogenetic tree was executed.
#
# The dated phylogenetic tree will be automatically created by the Nextflow
# pipeline and will be located in the results/add_duplicates directory. The rest
# of the files must be supplied manually, i.e. they must be present in the
# working directory.

library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to the project directory."
  ),
 make_option(
    c("-f", "--file"),
    type = "character",
    help = "Path to aci prediction results."
  ),
  make_option(
    c("-t", "--trees"),
    type = "character",
    help = "All dated trees collected in a single rds file."
  ),
  make_option(
    c("-s", "--sensitivity"),
    type = "character",
    help = "Path to phage sensitivity lab results."
  ),
  make_option(
    c("-C", "--fig3C"),
    type = "character",
    help = "Fig3C rds file."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "aci",
    file = "results/filter_crab/aci_collapse_geodate_crab.rds",
    trees = "data/all_dated_trees.rds",
    sensitivity = "labdata/aci_project/Ab_all_strains_phages_spotassay_PFU.tsv",
    fig3C = "results/plot_Fig3C_phylodist_phagedist/Fig3C_phylodist_phagedist.rds"
  )
}


library(devtools)
library(dplyr)
library(treeio)

load_all(args$project_dir)
load_all(paste0(args$project_dir, "/prophyl"))

# path to the phylogenetic tree
trees_path <- args$trees

# path to the file used for creating the phylogenetic tree. This file also
# contains metadata required by this script.
meta_path <- args$file

trees <- readRDS(trees_path)
meta_all <- readRDS(meta_path)

# Capitalise the first letters of country names
meta_all$country <- stringr::str_to_title(meta_all$country)
# Edit some country names manually
meta_all$country <- gsub("Bosnia_and_herzegovina", "Bosnia and Herzegovina", meta_all$country)
meta_all$country <- gsub("United_kingdom", "UK", meta_all$country)

# path to the file which contains assemblies that should be highlighted
labres_path <- args$sensitivity
labres_all <- read.csv(args$sensitivity, sep = "\t")

# filter to unique assemblies
labres_all <- labres_all[,c("Assembly", "MLST", "KL")]
labres_all <- dplyr::distinct(labres_all)
labres_all <- labres_all[complete.cases(labres_all),]

# import country colors
countries <- read.csv(
  paste0(args$project_dir, "/data/geographic_locations_in_study.tsv"),
  sep = "\t"
) %>%
  filter(variable == "country" & colour != "#000000")

country_colors_all <- data.frame(
  country = countries$term_pretty,
  color = countries$colour
)

city_colors_all <- data.frame()

serotrees <- list()

fig3C <- readRDS(args$fig3C)

fig3C_serotypes <- unique(fig3C$data$serotype)

tree_stats <- data.frame()

for (t in seq_along(trees)) {
  tree <- trees[[t]]
  ST <- names(trees)[t]
  
  if (length(ST) > 1) {
    stop("Script is designed to work with a single ST.")
  }
  
  meta <- meta_all[which(meta_all$mlst == ST),]
  
  labres <- labres_all[which(labres_all$MLST == ST),]
  
  # Throw a message if tested assemblies are missing from the tree
  
  # missing assemblies
  missing_assemblies <- labres$Assembly[which(
    labres$Assembly %in% meta$assembly == FALSE
  )]
  
  # print informative message
  if (length(missing_assemblies) > 0) {
    missing_assemblies_collapsed <- paste(
      missing_assemblies, collapse = ", "
    )
    msg <- paste0(
      "Some assemlies are missing from the tree: ",
      missing_assemblies_collapsed, "."
    )
    message(msg)
  }
  
  # drop tips which do not have related metadata
  if (all(tree$tip.label %in% meta$assembly) == FALSE) {
    labels <- tree$tip.label[which(tree$tip.label %in% meta$assembly == FALSE)]
    msg <- paste(labels, collapse = ", ")
    warning("One or more tip labels cannot be found in metadata table: ", msg)
    # drop tips that are not in the table
    tree <- ape::drop.tip(tree, tip = labels)
  }
  
  # only keep metadata for assemblies that are on the tree
  meta <- meta[which(meta$assembly %in% tree$tip.label),]
  
  k_types <- unique(labres$KL)
  
  # prepare a tree for each K type
  for (k in k_types) {
    
    sero <- paste0(ST, " - ", k)

    print(sero)
    
    
    # subset to those with given K type
    # subset to those from Europe
    # subset to those that are on the tree
    meta_small <- meta[which(
      meta$k_serotype == k &
        meta$continent == "europe" &
        meta$assembly %in% tree$tip.label
    ),]
    small_tree <- ape::keep.tip(tree, meta_small$assembly)
    small_tree_tbl <- treeio::as_tibble(small_tree)
    small_tree_tbl <- dplyr::left_join(
      small_tree_tbl,
      meta_small[,c(
        "assembly", "continent", "country", "city", "mlst",
        "k_serotype"
      )],
      by = c("label" = "assembly")
    )
    small_tree_tbl$tested <- FALSE
    small_tree_tbl$tested[which(small_tree_tbl$label %in% labres$Assembly)] <- TRUE
    s1 <- 0
    s2 <- 1
    while (s2 > s1) {
      s1 <- sum(small_tree_tbl$tested, na.rm = TRUE)
      parents <- unique(small_tree_tbl$parent[which(small_tree_tbl$tested == TRUE)])
      index <- which(small_tree_tbl$node %in% parents)
      small_tree_tbl$tested[small_tree_tbl$node %in% index] <- TRUE
      s2 <- sum(small_tree_tbl$tested, na.rm = TRUE)
    }
    # use single color for MLST
    mlst_colors <- data.frame(
      mlst = ST,
      color = "grey50"
    )
    # use single color for K type
    k_colors <- data.frame(
      k_serotype = k,
      color = "grey50"
    )
    # define colors for countries
    
    country_colors <- get_colors(
      df = small_tree_tbl,
      var = "country",
      colors = country_colors_all,
      verbose = TRUE
    )
    
    country_colors_all <- dplyr::bind_rows(
      country_colors_all,
      country_colors
    ) %>% distinct()
    
    # define colors for cities
    # only include cities which come from countries with multiple cities
    # citycount <- small_tree_tbl %>% 
    #   group_by(country) %>% 
    #   summarise(count = length(unique(city[which(!is.na(city))])))
    # selected_countries <- citycount$country[which(citycount$count > 1)]
    # if (length(selected_countries) == 0) {
    #   index <- numeric()
    # } else {
    #   index <- which(small_tree_tbl$country %in% selected_countries)
    # }
    # if (length(index) == 0) {
    # } else {
    #   selected_cities <- unique(small_tree_tbl$city[index])
    #   selected_cities <- selected_cities[which(!is.na(selected_cities))]
    #   selected_cities_colors <- qualpalr::qualpal(length(selected_cities), "pretty")$hex
    # }
    # city_colors <- data.frame(
    #   city = sort(selected_cities),
    #   color = selected_cities_colors
    # )
    
    city_colors <- get_colors(
      df = small_tree_tbl,
      var = "city",
      colors = city_colors_all,
      na.color = "#808080",
      verbose = TRUE
    )
    
    city_colors_all <- dplyr::bind_rows(
      city_colors_all,
      city_colors
    ) %>% distinct()
    
    # for ST492 use
    # - heatmap_width = 0.75
    # for all other STs, use
    # - heatmap_width = 2.5

    tree_stats <- dplyr::bind_rows(
      tree_stats,
      data.frame(
        serotype = sero,
        Ntips = ape::Ntip(small_tree),
        heatmap_width = dplyr::case_when(
          sero == "ST1 - KL1" ~ 6.5,
          sero == "ST1 - KL17" ~ 6.5,
          sero == "ST2 - KL2" ~ 5.5,
          sero == "ST2 - KL3" ~ 2.5,
          sero == "ST2 - KL9" ~ 5.6,
          sero == "ST2 - KL12" ~ 3.8,
          sero == "ST2 - KL32" ~ 3.5,
          sero == "ST2 - KL77" ~ 5.5,
          sero == "ST636 - KL40" ~ 4.2,
          sero == "ST492 - KL104" ~ 0.95,
          .default = 2.5
        )
      )
    )

    hmwidth <- tree_stats %>% filter(serotype == sero) %>% pull(heatmap_width)

    g <- plot_tree_fan(
      small_tree_tbl,
      linewidth = 0.1,
      highlight_var = "tested",
      open_angle = 15,
      heatmap_var = c("country", "city"),
      heatmap_colors = list(
        "country" = country_colors,
        "city" = city_colors
      ),
      legend_show = c("country", "tested"),
      heatmap_width = hmwidth,
      heatmap_colnames_font_size = 2,
      verbose = TRUE
    )
    
    if (sero %in% fig3C_serotypes) {
      serotrees[[sero]] <- g
      names(serotrees)[length(serotrees)] <- sero
    }
    
    # Export each tree as pdf

    if (!dir.exists("serotype_trees")) {
      dir.create("serotype_trees")
    }

    if (!sero %in% c(
      "ST2 - KL3",
      "ST636 - KL40",
      "ST492 - KL104",
      "ST2 - KL2",
      "ST1 - KL1",
      "ST2 - KL9",
      "ST2 - KL12",
      "ST1 - KL17",
      "ST2 - KL7",
      "ST2 - KL32",
      "ST2 - KL77"
    )) next()
    
    g <- g + ggtitle(paste0(ST, " - ", k, " Europe"))

    ggsave(
      filename = paste0(
        "serotype_trees/",
        ST, "_", k,
        "_Europe.pdf"
      ),
      plot = g,
      limitsize = FALSE
    )

    ggsave(
      filename = paste0(
        "serotype_trees/",
        ST, "_", k,
        "_Europe.png"
      ),
      plot = g,
      width = 8,
      height = 6
    )

    g <- g + geom_tiplab2(align = TRUE, offset = 10, size = 1)
    
    ggsave(
      filename = paste0(
        "serotype_trees/",
        ST, "_", k,
        "_Europe_tiplabs.pdf"
      ),
      plot = g,
      limitsize = FALSE
    )

  }
}

library(patchwork)

g <- wrap_plots(serotrees, ncol = 2) + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = list(names(serotrees))) & 
  theme(
    plot.tag = element_text(size = 5),
    plot.margin = unit(c(0,0,0,0), "cm"),
    legend.margin = margin(0,0,0,0),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.key.size = unit(2, "mm")
  )

ggsave(
  filename = "FigS6_serotype_trees.pdf",
  plot = g,
  width = 18,
  height = 13,
  units = "cm"
)

ggsave(
  filename = "FigS6_serotype_trees.png",
  plot = g,
  limitsize = FALSE,
  width = 18,
  height = 13,
  units = "cm"
)

write.table(
  tree_stats,
  file = "tree_stats.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
