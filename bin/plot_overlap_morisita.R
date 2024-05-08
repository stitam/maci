library(dplyr)
library(ggplot2)
library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-f", "--file_path"),
    type = "character",
    help = "Path to the input file, e.g. country_comparisons_collapse_geodate.tsv."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    file_path = "results/calc_serotype_freqs/geodate/country_comparisons_collapse_geodate.tsv"
  )
}

xtext_1 <- "within Europe"
xtext_2 <- "between continents"

data<- read.csv(args$file_path, sep = "\t")

## remove pairs belonging to the "same continent, other than Europe" category

data2 <- data %>% 
         filter(case_when
         (country1_continent=="europe" & country2_continent=="europe" ~ country1_continent == country2_continent,
           TRUE ~ country1_continent != country2_continent))

# categorization
tbl  <- data2 %>%
    mutate(category=case_when(
    country1_continent=="europe" & country2_continent=="europe" ~ xtext_1,
    country1_continent != country2_continent ~ xtext_2
  ))

# plot the mean prevalence of overlapping serotypes
  
g <- tbl  %>% 
  ggplot(mapping = aes(x=category, y= 100*(mean_overlap_prevalence)))+
  theme_bw()+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)))+
  #theme(panel.background = element_rect(fill = "white"))+
  scale_x_discrete(limits=c(xtext_1,xtext_2))+
  geom_boxplot(outlier.color="transparent")+
  geom_point(position = position_jitter(0.1, seed = 0)) +
  ggpubr::geom_signif(
    comparisons=list(
      c(xtext_1, xtext_2)),
    map_signif_level=TRUE,
    step_increase = 0.1,
    test="wilcox.test"
  ) +
  xlab("") +
  ylab("Mean prevalence of overlapping serotypes \nbetween country pairs")

ggsave(
  filename = "Mean_prevalence_overlapping_serotypes_b.pdf",
  plot = g,
  width = 8,
  height = 6
)
ggsave(
  filename = "Mean_prevalence_overlapping_serotypes_b.png",
  plot = g,
  width = 8,
  height = 6
)
saveRDS(g, "Mean_prevalence_overlapping_serotypes_b.rds")

# plot Morisita index 

g <- tbl  %>% 
  ggplot(mapping = aes(x=category, y= morisita))+
  theme_classic()+
  theme(
    axis.text = element_text(size = 5),
    axis.title.y = element_text(
      size = 5,
      margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.line=element_line(color="black", linewidth=0.1),
    panel.grid.major.y=element_line(color="gray", linetype = 3, linewidth = 0.1),
    plot.margin = unit(c(0,0,0,0), "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
  ) +
  #theme(panel.background = element_rect(fill = "white"))+
  scale_x_discrete(limits=c(xtext_1,xtext_2))+
  geom_boxplot(outlier.color="transparent", lwd = 0.1)+
  geom_point(position = position_jitter(0.1, seed = 0), size = 0.5) +
  ggpubr::geom_signif(
    comparisons=list(
      c(xtext_1, xtext_2)),
    map_signif_level=function(p) {
      g$morisita_p <<- p
      if (p < 0.001) {
        return("***")
      } else if (p < 0.01) {
        return("**")
      } else if  (p < 0.05) {
        return("*")
      } else {
        return("NS")
      }
    },
    step_increase = 0.1,
    test="wilcox.test",
    size = 0.1,
    textsize = 2
  ) +
  xlab("") +
  ylab("Population composition dissimilarity \nbetween country pairs \n(Morisita index)") +
  ylim(0, 1.2)

ggsave(
  filename = "Fig1E_Morisita_countries.pdf",
  plot = g,
  width = 8,
  height = 6
)
ggsave(
  filename = "Fig1E_Morisita_countries.png",
  plot = g,
  width = 8,
  height = 6
)
saveRDS(g, "Fig1E_Morisita_countries.rds")

export_data <- g$data %>%
  select(
    country1,
    country2,
    country1_continent,
    country2_continent,
    morisita,
    category
  )

write.table(
  export_data,
  file = "Fig1E_Morisita_countries.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)