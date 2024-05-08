library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-r", "--region_file"),
    type = "character",
    help = "Path to a tbl prepared for serotop by regions, e.g. top_serotypes_region23_collapse_geodate.tsv."
  ),
  make_option(
    c("-s", "--serotype_file"),
    type = "character",
    help = "Path to a tbl which lists global or prevalent serotypes."
  ),
  make_option(
    c("-c", "--country_file"),
    type = "character",
    help = "Path to a tbl prepared for serotop by countries, e.g. top_serotypes_country_collapse_geodate.tsv."
  ),
  make_option(
   c("-S", "--shortlist"),
   type = "character",
   help = "Path to shortlisted regions and countries."
  ),
  make_option(
   c("-m", "--metadata"),
   type = "character",
   help = "Path to an rds file containing metadata."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    region_file = "results/calc_serotype_freqs/geodate/TableS2A_top_serotypes_region23_collapse_none.tsv",
    serotype_file = "results/calc_serotype_freqs/geodate/global_or_prevalent_serotypes_region23_collapse_none.tsv",
    country_file = "results/calc_serotype_freqs/geodate/top_serotypes_country_collapse_none.tsv",
    shortlist = "aci/data/geographic_locations_in_study.tsv",
    metadata = "data/typing_summary_tables/aci_study.rds"
  )
}

library(tidyverse)
library(ComplexHeatmap)

strategy <- args$region_file %>% basename() %>% strsplit(., "_") %>% unlist()
strategy <- gsub(".tsv", "", strategy[6])

year_range <- readRDS(args$metadata)

minyear <- year_range$minyear
maxyear <- year_range$maxyear

# GEODATE ###########
# REGIONS ######
regions_serotypes_geodate <- read_delim(args$region_file)
regions_serotypes_geodate$serotype <- gsub("_"," ", regions_serotypes_geodate$serotype)
regions_serotypes_geodate$region23 <- gsub("_"," ", regions_serotypes_geodate$region23)
regions_serotypes_geodate$region23 <- tools::toTitleCase(regions_serotypes_geodate$region23)

index_remaining <- which(regions_serotypes_geodate$ratio_cumsum > 0.9)
remaining_count <- length(unique(regions_serotypes_geodate$serotype[index_remaining]))

regions_serotypes_geodate_hm <- regions_serotypes_geodate  %>%  
  mutate(
    regions_serotypes_geodate,
    serotype2=if_else(
      ratio_cumsum>0.9, 
      paste0("Other (", remaining_count, ")"), serotype)) %>% 
  group_by(region23, serotype2) %>% 
  summarise(percent=sum(ratio*100))

regions_serotypes_geodate_wide <- regions_serotypes_geodate_hm %>% pivot_wider(names_from=serotype2, values_from = percent)
regions_serotypes_geodate_wide <- subset(regions_serotypes_geodate_wide, select = c(1, 3:ncol(regions_serotypes_geodate_wide), 2))
m <- as.matrix(regions_serotypes_geodate_wide[,-1])
rownames(m) <- regions_serotypes_geodate_wide$region23
cn = colnames(regions_serotypes_geodate_wide[-1])


# Bray-Curtis
regions_serotypes_geodate_bray <- regions_serotypes_geodate  %>%  
  mutate(regions_serotypes_geodate, serotype2=serotype) %>% 
  group_by(region23, serotype2) %>% 
  summarise(count=sum(count))

regions_serotypes_geodate_bray_wide <- pivot_wider(regions_serotypes_geodate_bray, names_from=serotype2, values_from=count)
diversity_matrix <- as.matrix(regions_serotypes_geodate_bray_wide[,-1])
rownames(diversity_matrix) <- regions_serotypes_geodate_bray_wide$region23
diversity_matrix[is.na(diversity_matrix)] <- 0
beta_diversity <- vegan::vegdist(diversity_matrix, "bray")

# dendogram
my_clustering <- hclust(dist(beta_diversity))


##### EDIT #####
# Instead of defining global and prevalent clones here, define them in an
# upstream process and pass them in as an argument.

# coloring global and prevalent clones
# global_clones <- regions_serotypes_geodate %>% filter(ratio>0.02) %>% 
#   group_by(serotype) %>% 
#   summarise(n_region=length(unique(region23))) %>% 
#   filter(n_region>=3) %>% subset(select=1) %>% filter(serotype!="ST25 KL14") %>% mutate(color="#DC143C")

# prevalent_clones <- regions_serotypes_geodate %>% filter(ratio>0.05) %>% summarise(serotype=unique(serotype)) %>% mutate(color="orange3")

# prevalent_clones <- anti_join(prevalent_clones, global_clones, by="serotype")
# global_prevalent <- rbind(global_clones, prevalent_clones)

global_prevalent <- read.csv(args$serotype_file, sep = "\t")
global_prevalent$serotype <- gsub("_", " ", global_prevalent$serotype)

################

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

# TODO make this OS independent
#library(extrafont)
#font_import(
#  path = "/usr/share/fonts/truetype/msttcorefonts", 
#  prompt = FALSE
#)

index_prevalent <- which(sero_color$serotype %in% global_prevalent$serotype)
index_prevalent <- c(index_prevalent, nrow(sero_color))

xtext_linecolors <- rep(
  c("#1F1E53", "#925637"),
  times = ceiling(length(index_prevalent)/2)
)
xtext_linecolors <- xtext_linecolors[1:length(index_prevalent)]

heatmap_region <- Heatmap(
  m,
  cluster_rows = my_clustering,
  cluster_columns = F,
  row_dend_width = unit(0.4, "cm"),
  row_dend_gp = gpar(lwd = 0.1),
  row_names_gp = gpar(fontsize = 5),
  col=colors,
  heatmap_legend_param = list(
    title="Percent",
    title_gp = gpar(fontsize = 5),
    at = gbreaks,
    labels_gp = gpar(fontsize = 5),
    legend_height = unit(1, "cm"),
    legend_width = unit(2, "mm"),
    color_bar = "continuouss"
  ),
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  na_col="white",
  show_column_names = FALSE,
  row_names_side = "left",
  bottom_annotation = HeatmapAnnotation(
    text = anno_mark(
      at = index_prevalent,
      labels = sero_color$serotype[index_prevalent],
      which = "column",
      side = "bottom",
      lines_gp = gpar(lwd = 0.5, col = xtext_linecolors ),
      labels_gp = gpar(fontsize = 5, col=sero_color$color[index_prevalent]),
      labels_rot = 60,
      padding = unit(5, "mm")
    )
  )
) 

heatmap_region

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

ggsave(
    g1,
    filename = paste0("Fig1D_heatmap_region23_collapse_", strategy, ".pdf"),
    height=5,
    width=20
)
ggsave(
    g1,
    filename = paste0("Fig1D_heatmap_region23_collapse_", strategy, ".png"),
    height=5,
    width=20
)
saveRDS(
    g1,
    file = paste0("Fig1D_heatmap_region23_collapse_", strategy, ".rds")
)

# vertical
m_t <- t(m)
rownames(m_t) <- colnames(regions_serotypes_geodate_wide[-1])
cn_m_t = regions_serotypes_geodate_wide$region23

# breaks for plot legend
gbreaks <- c(2,5)
maxm <- 10 * floor(max(m_t, na.rm = TRUE) / 10)
if (maxm >= 0) {
  gbreaks <- c(gbreaks, seq(from = 10, to = maxm, by = 10))
}

heatmap_region_vertical <- Heatmap(m_t, cluster_columns = my_clustering, cluster_rows = F,
                           col=colors,
                           heatmap_legend_param = list(title="Percent", at = gbreaks, labels = gbreaks),
                           rect_gp = gpar(col = "grey80", lwd = 0.5), na_col="white", show_column_names = FALSE, row_names_side = "left",
                           bottom_annotation = HeatmapAnnotation(
                            text = anno_text(cn_m_t, rot=30, location = unit(1, "npc"), just = "right",gp = gpar(fontsize = 8))),
                           row_names_gp = gpar(fontsize = 8, col=sero_color$color))

g2 <- ggplotify::as.ggplot(heatmap_region_vertical)
ggsave(
    g2,
    filename = paste0(
        "Fig1D_heatmap_region23_vertical_collapse_", strategy, ".pdf"),
    height=15,
    width=5
)
ggsave(
    g2,
    filename = paste0(
        "Fig1D_heatmap_region23_vertical_collapse_", strategy, ".png"),
    height=15,
    width=5
)
saveRDS(
    g2,
    file = paste0("Fig1D_heatmap_region23_vertical_collapse_", strategy, ".rds")
)

# COUNTRIES ######
country_serotypes_geodate <- read_delim(args$country_file)
country_serotypes_geodate$serotype <- gsub("_"," ", country_serotypes_geodate$serotype)
country_serotypes_geodate$country <- gsub("_"," ", country_serotypes_geodate$country)
country_serotypes_geodate$country <- tools::toTitleCase(country_serotypes_geodate$country)


# these countries do not match fully the shortlisted european countries
# TODO get this list from the shortlist file.
shortlist <- read.csv(args$shortlist, sep = "\t") %>% 
  dplyr::filter(variable == "country" &  shortlisted == TRUE)

index_europe <- which(
  countrycode::countrycode(
    shortlist$term,
    origin = "country.name",
    destination = "continent") == "Europe"
)

countries <- shortlist$term_pretty[index_europe]

country_serotypes_geodate_hm <- country_serotypes_geodate  %>%  
  filter(country_serotypes_geodate$country %in% countries) %>% 
  group_by(country, serotype) %>% 
  summarise(percent=sum(ratio*100))

country_serotypes_geodate_wide <- country_serotypes_geodate_hm %>% pivot_wider(names_from=serotype, values_from = percent)
m2 <- as.matrix(country_serotypes_geodate_wide[, -1])
rownames(m2) <- country_serotypes_geodate_wide$country
cn2 = colnames(country_serotypes_geodate_wide[-1])


# Bray-Curtis
country_serotypes_geodate_bray <- country_serotypes_geodate  %>%  
  filter(country_serotypes_geodate$country %in% countries) %>% 
  group_by(country, serotype) %>% 
  summarise(count=sum(count))

country_serotypes_geodate_bray_wide <- pivot_wider(country_serotypes_geodate_bray, names_from=serotype, values_from=count)
diversity_matrix2 <- as.matrix(country_serotypes_geodate_bray_wide[,-1])
rownames(diversity_matrix2) <- country_serotypes_geodate_bray_wide$country
diversity_matrix2[is.na(diversity_matrix2)] <- 0
beta_diversity2 <- vegan::vegdist(diversity_matrix2, "bray")

# dendogram
my_clustering2 <- hclust(dist(beta_diversity2))

rest2 <- data.frame(serotype=colnames(country_serotypes_geodate_wide[-1])) %>% mutate(color="black")
rest2 <- anti_join(rest2,global_prevalent,  by="serotype")
sero_color2 <- rbind(global_prevalent, rest2)
sero_color2 <- sero_color2[match(colnames(country_serotypes_geodate_wide[-1]), sero_color2$serotype),]

# breaks for plot legend
gbreaks <- c(0,2,5)
maxm <- 10 * floor(max(m2, na.rm = TRUE) / 10)
if (maxm >= 0) {
  gbreaks <- c(gbreaks, seq(from = 10, to = maxm, by = 10))
}

heatmap_country <- Heatmap(m2, cluster_rows = my_clustering2, cluster_columns = F,
        col=colors,
        heatmap_legend_param = list(title="Percent", at = gbreaks, labels = gbreaks),
        rect_gp = gpar(col = "grey80", lwd = 0.5), na_col="white", show_column_names = FALSE, row_names_side = "left",
        bottom_annotation = HeatmapAnnotation(
          text = anno_text(cn2, rot = 60, location = unit(1, "npc"), just = "right",gp = gpar(fontsize = 8, col=sero_color2$color))))

g3 <- ggplotify::as.ggplot(heatmap_country) + 
  geom_label(
    aes(x = 0.85, y = 0.85, label = paste0(minyear, "-", maxyear)),
    fill = "#FFFFFF"
  )

g3$labels$label <- paste0(minyear, "-", maxyear)
g3$layers[[3]]$constructor[[2]]$label <- paste0(minyear, "-", maxyear)
g3$layers[[3]]$mapping$label <- paste0(minyear, "-", maxyear)

ggsave(
    g3,
    filename = paste0(
        "FigS2C_heatmap_country_collapse_", strategy, ".pdf"),
    height=4,
    width=10
)
ggsave(
    g3,
    filename = paste0(
        "FigS2C_heatmap_country_collapse_", strategy, ".png"),
    height=4,
    width=10
)
saveRDS(
    g3,
    file = paste0("FigS2C_heatmap_country_collapse_", strategy, ".rds")
)
