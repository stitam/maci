library(dplyr)
library(ggplot2)
rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
aci_path <- args[1]
pp_path <- args[2]

# import data set
aci <- readRDS(aci_path)

# add poppunk clusters
pp <- read.csv(pp_path)
pp <- dplyr::rename(pp, assembly = Taxon, pp = Cluster)
pp$assembly <- gsub("_", ".", pp$assembly)
pp$assembly <- gsub("GCA\\.", "GCA_", pp$assembly)
pp$assembly <- gsub("GCF\\.", "GCF_", pp$assembly)
aci <- dplyr::left_join(aci, pp, by = "assembly")

# filter to genomes that have predicted poppunk clusters
aci <- aci[which(!is.na(aci$pp)), ]

# calculate summary statistics
stats <- aci %>% group_by(pp) %>% summarise(
  count = n(),
  mlst = length(unique(mlst)),
  k_serotype = length(unique(k_serotype)),
  serotype_count  = length(unique(serotype)),
  serotype_list = paste(unique(serotype), collapse = ","),
  country = length(unique(country)),
  region23 = length(unique(region23)),
  continent = length(unique(continent)),
  collection_day_min = min(collection_day, na.rm = TRUE),
  collection_day_max = max(collection_day, na.rm = TRUE)
)

# format variables
stats$k_serotype <- factor(stats$k_serotype, levels = sort(unique(stats$k_serotype)))
stats$country <- factor(stats$country, levels = sort(unique(stats$country)))
stats$country <- factor(stats$country, levels = sort(unique(stats$country)))
stats$region23 <- factor(stats$region23, levels = sort(unique(stats$region23)))
stats$continent <- factor(stats$continent, levels = sort(unique(stats$continent)))

# define new variable
stats$yeardiff <- ifelse(
  is.infinite(stats$collection_day_max) | is.infinite(stats$collection_day_min),
  NA,
  stats$collection_day_max-stats$collection_day_min
)
stats$yeardiff <- signif(stats$yeardiff/365.25, 2)

saveRDS(stats, "stats_by_pp.rds")

write.table(
    stats,
    file = "stats_by_pp.tsv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

g <- ggplot(stats, aes(count)) + geom_histogram() + scale_x_log10()

ggsave("poppunk_cluster_sizes.png")

# only consider clusters that have at least 5 members
statsA <- stats[which(stats$count >= 5),]

g <- ggplot(stats, aes(k_serotype)) + 
  geom_bar(width = 0.5) + 
  scale_y_log10() +
  xlab("Number of K serotypes within a poppunk cluster") +
  ylab("Number of poppunk clusters")

ggsave("number_of_K_serotypes_within_a_poppunk_cluster.png")

g <- ggplot(stats, aes(country)) + 
  geom_bar(width = 0.5) + 
  scale_y_log10() +
  xlab("Number of countries within a poppunk cluster") +
  ylab("Number of poppunk clusters")

ggsave("number_of_countries_within_a_poppunk_cluster.png")

aci_summary <- list(
    "number_of_isolates" = sum(stats$count),
    "number_of_pp_clusers" = nrow(stats),
    "number_of_pp_clusters_count_min_2" = length(which(stats$count > 1))
)

saveRDS(aci_summary, "aci_summary.rds")