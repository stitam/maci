library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-A", "--figA_pdf"),
    type = "character",
    help = ""
  ),
  make_option(
    c("-a", "--figA_png"),
    type = "character",
    help = ""
  ),
  make_option(
    c("-B", "--figB"),
    type = "character",
    help = ""
  ),
  make_option(
    c("-C", "--figC"),
    type = "character",
    help = ""
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    figA_pdf = "./results_redacted/plot_sensitivity_heatmap/heatmap_phages_PFU_composed_country_order_noLfw_cities_old_coloring_pdf.rds",
    figA_png = "./results_redacted/plot_sensitivity_heatmap/heatmap_phages_PFU_composed_country_order_noLfw_cities_old_coloring_png.rds",
    figB = "./results_redacted/plot_Fig3B_sensitivity_dissimilarity/Fig3B_sensitivity_dissimilarity.rds",
    figC = "./results_redacted/plot_phylodist_phagedist/phylodist_phagedist.rds"
  )
}

library(ggplot2)
library(patchwork)

figA_pdf <- readRDS(args$figA_pdf) |> wrap_plots()
figA_png <- readRDS(args$figA_png) |> wrap_plots()
figB <- readRDS(args$figB)
figC <- readRDS(args$figC)

prep_fig3 <- function(figA, figB, figC) {
  figA <- figA + theme(
    axis.text = element_text(size = 5),
    strip.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5)
  )

  figB <- figB + theme(
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    axis.ticks = element_line(size = 0.1),
    panel.grid.major.y=element_line(color="gray", linetype = 3, linewidth = 0.1)
  )

  figC <- figC + theme_classic() + theme(
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    axis.line = element_line(color="black", linewidth=0.1),
    axis.ticks = element_line(size = 0.1),
    strip.text = element_text(size = 5),
    strip.background = element_rect(fill = "transparent", color = "black", linewidth = 0.1),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.margin = unit(c(0,0,0,0), "cm"),
    legend.box.margin = unit(c(0,0,0,0), "cm"),
    legend.position = "bottom",
    panel.background = element_rect(fill = "transparent", color = "black", linewidth = 0.1),
    panel.grid.major.y=element_line(color="gray", linetype = 3, linewidth = 0.1)
  )

  fig <- figA + figB + figC + 
    plot_layout(
      design = c(
        # row start, col start, row end, col end
        area(1, 1, 50, 100),
        area(51, 1, 100, 50),
        area(51, 51, 100, 100)
      )
    ) +
    plot_annotation(tag_levels = list(c("A", "", "", "B", "C"))) & 
    theme(
      plot.tag = element_text(size = 5)
    )
  return(fig)
}

fig_pdf <- prep_fig3(figA_pdf, figB, figC)
fig_png <- prep_fig3(figA_png, figB, figC)

ggsave(
  filename = "Fig4.png",
  plot = fig_png,
  width = 18,
  height = 12,
  units = "cm"
)

ggsave(
  filename = "Fig4.pdf",
  plot = fig_pdf,
  width = 18,
  height = 12,
  units = "cm"
)
