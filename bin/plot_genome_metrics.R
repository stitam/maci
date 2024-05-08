# This script is used to prepare the plots for the genome metrics of isolates.
# Genome metrics are compared between genomes sequenced in our lab and genomes
# downloaded from NCBI.
#
# TODO: the script filters genomes by QC and CRAB status. Check if the processes
# "collapse_outbreaks" and "filter_crab" could be swapped. This would allow this
# script to be run after "filter_crab" and would yield a much cleaner workflow.

library(optparse)
library(ggpubr)
library(qualpalr)
library(scales)
library(tidyverse)

rm(list = ls())

args_list <- list(
  make_option(
    c("-f", "--file"),
    type = "character",
    help = "Path to aci prediction results."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    file = "results/filter_assemblies/aci_filtered.rds"
  )
}

aci_all <- readRDS(args$file)

aci_all$source <- aci_all$assembly_source

# THIS IS WHERE GENOMES ARE FILTERED BY QC AND CRAB STATUS
aci_qc_true_crab <- aci_all %>% filter(qc_pass=="TRUE" & crab=="TRUE") 

set.seed(555)

nsource <- length(unique(aci_all$source))
if (nsource == 1) {
  colors <- data.frame(
    "source" = unique(aci_all$source),
    "color" = "black"
  )
} else if (nsource == 2) {
  colors <- data.frame(
    "source" = unique(aci_all$source),
    "color" = c("brown4", "deepskyblue4")
  )
} else {
  colors <- data.frame(
    "source" = sort(unique(aci_all$source)),
    "color" = qualpalr::qualpal(length(unique(aci_all$source)), "pretty")$hex
  )
}

color_vector <- colors$color
names(color_vector) <- colors$source

# contig number
contig_number <- ggplot(aci_qc_true_crab)+
  stat_boxplot(aes(source, contig_count), geom='errorbar', linetype=1, width=0.2)+
  geom_boxplot(aes(source, contig_count, fill=source), outlier.colour = "black",
               outlier.fill = "white")+
  ylab("Number of contigs")+
  xlab("Source of genome sequence")+
  scale_fill_manual(values = color_vector)+
  theme_bw()+
  theme(panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        legend.position = "none")


# N50
n50 <- ggplot(aci_qc_true_crab)+
  stat_boxplot(aes(source, N50), geom='errorbar', linetype=1, width=0.2)+
  geom_boxplot(aes(source, N50, fill=source), outlier.colour = "black",
               outlier.fill = "white")+
  ylab("N50 (bp)")+
  xlab("Source of genome sequence")+
  scale_fill_manual(values = color_vector)+
  theme_bw()+
  theme(panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        legend.position = "none")+
  scale_y_log10(labels=trans_format("log10", math_format(10^.x)),
                breaks=c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7))

# N95
n95 <- ggplot(aci_qc_true_crab)+
  stat_boxplot(aes(source, N95), geom='errorbar', linetype=1, width=0.2)+
  geom_boxplot(aes(source, N95, fill=source), outlier.colour = "black",
               outlier.fill = "white")+
  ylab("N95 (bp)")+
  xlab("Source of genome sequence")+
  scale_fill_manual(values = color_vector)+
  theme_bw()+
  theme(panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        legend.position = "none")+
  scale_y_log10(labels=trans_format("log10", math_format(10^.x)),
                breaks=c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7))


# genome size
genome_size <- ggplot(aci_qc_true_crab)+
  stat_boxplot(aes(source, genome_size/10^6), geom='errorbar', linetype=1, width=0.2)+
  geom_boxplot(aes(source, genome_size/10^6, fill=source), outlier.colour = "black",
               outlier.fill = "white")+
  ylab("Genome size (Mbp)")+
  xlab("Source of genome sequence")+
  scale_fill_manual(values = color_vector)+
  theme_bw()+
  theme(panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        legend.position = "none")


# gc content
gc_content <- ggplot(aci_qc_true_crab)+
  stat_boxplot(aes(source, gc_content*100), geom='errorbar', linetype=1, width=0.2)+
  geom_boxplot(aes(source, gc_content*100, fill=source), outlier.colour = "black",
               outlier.fill = "white")+
  ylab("G+C conent (%)")+
  xlab("Source of genome sequence")+
  scale_fill_manual(values = color_vector)+
  theme_bw()+
  theme(panel.grid.major = element_line(color="white"),
        panel.grid.minor = element_line(color="white"),
        legend.position = "none")

genome_metrics <- ggarrange(contig_number, n50, n95, genome_size, gc_content, ncol=3, nrow=2)

ggsave(
  "FigS1_genome_metrics.pdf",
  width = 8,
  height = 6
)

ggsave(
  "FigS1_genome_metrics.png",
  width = 8,
  height = 6
)
