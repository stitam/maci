library(ggplot2)
library(ggpubr)
library(optparse)
library(patchwork)
rm(list = ls())

args_list <- list(
  make_option(
    c("-A", "--fig1A"),
    type = "character",
    help = ""
  ),
  make_option(
    c("-B", "--fig1B"),
    type = "character",
    help = ""
  ),
  make_option(
    c("-C", "--fig1C"),
    type = "character",
    help = ""
  ),
  make_option(
    c("-D", "--fig1D"),
    type = "character",
    help = ""
  ),
  make_option(
    c("-E", "--fig1E"),
    type = "character",
    help = ""
  ),
  make_option(
    c("-r", "--regions"),
    type = "character",
    help = "Path to a table with region names, coordinates, colors."
  )
  
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    fig1A = "Fig1A_world_map.rds",
    fig1B = "Fig1B_crab_over_time.rds",
    fig1C = "Fig1C_serotop_region23_collapse_geodate.rds",
    fig1D = "Fig1D_heatmap_region23_collapse_geodate.rds",
    fig1E = "Fig1E_Morisita_countries.rds",
    regions = "geographic_locations_in_study.tsv"
  )
}

fig1A <- readRDS(args$fig1A)
fig1B <- readRDS(args$fig1B)
fig1C <- readRDS(args$fig1C)
fig1D <- readRDS(args$fig1D)
fig1E <- readRDS(args$fig1E)

fig1 <- fig1A + fig1B + fig1C + guide_area() + fig1D + fig1E +
  plot_layout(
    design = c(
      # row start, col start, row end, col end
      area(1, 1, 94, 35),
      area(1, 36, 90, 51),
      area(1, 52, 90, 70),
      area(95, 1, 110, 70),
      area(111, 1, 230, 52),
      area(111, 54, 225, 70)
    ),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") & 
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 5)
  )

ggsave(
  filename = "Fig1.png",
  plot = fig1,
  width = 18,
  height = 10,
  units = "cm"
)

ggsave(
  filename = "Fig1.pdf",
  plot = fig1,
  width = 18,
  height = 10,
  units = "cm"
)
