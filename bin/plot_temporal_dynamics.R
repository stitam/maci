library(extrafont)
library(ggplot2)
library(ggpubr)
library(optparse)
library(patchwork)
rm(list = ls())

# required for exporting proper fonts with pdf
#extrafont::loadfonts()
#extrafont::choose_font("Times New Roman")


args_list <- list(
  make_option(
    c("-A", "--figA"),
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
  ),
  make_option(
    c("-D", "--figD"),
    type = "character",
    help = ""
  ),
  make_option(
    c("-E", "--figE"),
    type = "character",
    help = ""
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    figA = "results_redacted/plot_sero_over_time/geodate/sero_over_time_ds_geodate.rds",
    figB = "results_redacted/plot_morisita_histograms/Morisita_histograms_geodate_crab.rds",
    figC = "results_redacted/plot_global_tree/geodate/global_tree.rds",
    figD = "data/rr_global_no_focus/relative_risks_type2.rds",
    figE = "data/rr_regional/relative_risks_type3.rds"
  )
}

figA <- readRDS(args$figA)

figA <- figA +
  guides(
    fill=guide_legend(title = "", nrow = 4)
  )

figB <- readRDS(args$figB)+
  theme(
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    strip.text = element_text(size = 5)
  ) #+ 
  #geom_text(aes(
  #  x = over_time_x,
  #  y = over_time_y,
  #  label = "2009-2015\nvs\n2016-2022",
  #), col = "red", size = 5/.pt, hjust = "left", vjust = "top",lineheight = .7) 

figC <- readRDS(args$figC)

figC <- figC +
  theme(
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.key.size = unit(2, "mm"),
    legend.margin = margin(0,0,0,0),
    plot.margin = unit(c(0,0,0,0), "cm")
  )

figD <- readRDS(args$figD) + theme(
  plot.margin = unit(c(0,0,0,0), "cm"),
  axis.text.x = element_text(angle = 0, hjust = 0.5)
)

figE <- readRDS(args$figE) + theme(
  plot.margin = unit(c(0,0,0,0), "cm"),
  axis.text.x = element_text(angle = 0, hjust = 0.5)
) + ylab("")

fig_pdf <- figA + figB + figC + figD + figE + guide_area() +
  plot_layout(
    design = c(
      # row start, col start, row end, col end
      area(1, 1, 84, 25),
      area(1, 28, 41, 58),
      area(1, 59, 84, 100),
      area(42, 28, 84, 45),
      area(42, 46, 84, 58),
      area(85, 1, 100, 100)
    ),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") & 
  theme(
    legend.position = "bottom",
    #legend.margin = unit(c(0,0,0,0), "cm"),
    plot.tag = element_text(size = 5)
  )

ggsave(
  filename = "Fig3.pdf",
  plot = fig_pdf,
  width = 18,
  height = 10,
  units = "cm"
)

fig_png <- figA + figB + figC + figD + figE + guide_area() +
  plot_layout(
    design = c(
      # row start, col start, row end, col end
      area(1, 1, 84, 25),
      area(1, 28, 41, 60),
      area(1, 61, 84, 100),
      area(42, 28, 84, 45),
      area(42, 46, 84, 60),
      area(85, 1, 100, 100)
    ),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") & 
  theme(
    legend.position = "bottom",
    #legend.margin = unit(c(0,0,0,0), "cm"),
    plot.tag = element_text(size = 5)
  )

ggsave(
  filename = "Fig3.png",
  plot = fig_png,
  width = 18,
  height = 10,
  units = "cm"
)
