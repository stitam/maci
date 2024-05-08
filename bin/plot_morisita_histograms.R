library(dplyr)
library(ggplot2)
library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-a", "--all_counts"),
    type = "character",
    help = "Path to file that includes serotype counts by country."
  ),
  make_option(
    c("-y", "--by_year"),
    type = "character",
    help = "Path to file that includes serotype counts by country and year."
  ),
  make_option(
    c("-c", "--comparisons"),
    type = "character",
    help= "Path to file that includes serotype prevalence comparisons between countries."
  ),
  make_option(
    c("-m", "--minyear"),
    type = "integer",
    help = "Lowest year to include."
  ),
  make_option(
    c("-r", "--minyear_recent"),
    type = "integer",
    help = "Lowest year to include for recent time period."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    all_counts = "results/calc_serotype_freqs/geodate/TableS2A_top_serotypes_country_collapse_geodate.tsv",
    by_year = "results/calc_serotype_freqs/geodate/serotypes_country_year_collapse_geodate.tsv",
    comparisons = "results/calc_serotype_freqs/geodate/country_comparisons_collapse_geodate.tsv",
    minyear = 2009,
    minyear_recent = 2016
  )
}

collapse_strategy <- args$by_year %>% basename() %>% strsplit(split = "_") %>% unlist()
collapse_strategy <- collapse_strategy[5]
collapse_strategy <- gsub("\\.tsv", "", collapse_strategy)

# define countries to select
countries <- c("china", "germany", "usa")

# import data
aci <- read.csv(args$by_year, sep = ",")

# keep selected countries
aci <- aci %>% filter(country %in% countries)

# construct bin names for recent time period and the period before
old_bin_name <- paste0("(", min(aci$collection_year)-1, ",", as.numeric(args$minyear)-1, "]")
non_recent_bin_name <- paste0("(", as.numeric(args$minyear)-1, ",", as.numeric(args$minyear_recent)-1, "]")
recent_bin_name <- paste0("(", as.numeric(args$minyear_recent)-1 , ",", max(aci$collection_year), "]")

# create new variable collection_year_bin
aci$collection_year_bin <- cut(
  aci$collection_year,
  breaks = c(
    0,
    as.numeric(args$minyear)-1,
    as.numeric(args$minyear_recent)-1,
    max(aci$collection_year)
  ),
  labels <- c(
    old_bin_name,
    non_recent_bin_name,
    recent_bin_name
  )
)

# filter to recent samples and to samples from the time period before
aci <- aci %>% filter(collection_year_bin %in% c(non_recent_bin_name, recent_bin_name))

# collapse by country, collection_year_bin, serotype
aci_collapsed <- aci %>% 
  group_by(country, collection_year_bin, serotype) %>% 
  summarise(count = sum(count))

# create composte id from country and collection_year_bin
aci_collapsed$country_year <- paste(
    aci_collapsed$country,
    aci_collapsed$collection_year_bin, 
    sep = "_"
)

diversity_matrix <- tidyr::pivot_wider(
  aci_collapsed,
  id_cols = country_year,
  names_from = serotype, 
  values_from = count
) %>%
  tibble::column_to_rownames("country_year") %>%
  replace(is.na(.), 0)

yearbins <- as.character(unique(aci_collapsed$collection_year_bin))

morisita <- data.frame(
  country = countries,
  yearbin_1 = yearbins[1],
  yearbin_2 = yearbins[2],
  morisita = NA
)
for (i in 1:nrow(morisita)) {
  rn_1 <- paste(morisita$country[i], morisita$yearbin_1[i], sep = "_")
  index_1 <- which(rownames(diversity_matrix) == rn_1)
  rn_2 <- paste(morisita$country[i], morisita$yearbin_2[i], sep = "_")
  index_2 <- which(rownames(diversity_matrix) == rn_2)
  morisita$morisita[i] <- abdiv::morisita(
      diversity_matrix[index_1, ],
      diversity_matrix[index_2, ]
  ) %>% signif(., digits = 4)
}

# export morisita index

write.table(
  morisita,
  file = paste0(
    "morisita_country_year_collapse_",
    collapse_strategy,
    "_crab.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# prepare histograms

comparisons <- read.csv(args$comparisons, sep = "\t")

all_countries <- sort(unique(c(comparisons$country1, comparisons$country2)))

df <- data.frame(
  country1 = rep(countries, each = length(all_countries)),
  country2 = rep(all_countries, times = length(countries))
)
df <- df[which(df$country1 != df$country2),]
df$morisita <- NA
df$over_time <- NA

for (i in 1:nrow(df)) {
  index <- which(comparisons$country1 == df$country1[i] & 
                   comparisons$country2 ==  df$country2[i])
  if (length(index) == 0) {
    index <- which(comparisons$country2 == df$country1[i] & 
                     comparisons$country1 ==  df$country2[i])
  }
  df$morisita[i] <- comparisons$morisita[index]
  df$over_time[i] <- morisita$morisita[which(morisita$country == df$country1[i])]
}

df$country1 <- gsub("china", "China", df$country1)
df$country1 <- gsub("germany", "Germany", df$country1)
df$country1 <- gsub("usa", "USA", df$country1)

df$over_time_x <- sapply(df$country1, function(x) {
  if (x == "China") return(0)
  if (x == "Germany") return(0)
  if (x == "USA") return(0)
})

df$over_time_y <- 4

g <- ggplot(df, aes(morisita)) + 
  geom_histogram(col = "black", binwidth = 0.1, size = 0.1) + 
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  facet_wrap(country1 ~ ., ncol = 1) +
  geom_vline(aes(xintercept = over_time), col = "red", linewidth = 0.1) +
  xlab("Population composition dissimilarity\n(Morisita index)") +
  ylab("Number of country pairs") +
  theme_minimal() +
  theme(
    panel.background = element_rect(colour = "#000000", linewidth = 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(colour = "#000000", linewidth = 0.1),
    axis.ticks.length = unit(0.05, "cm"),
    strip.text = element_text(margin = margin(0,0,0,0, "cm"))
  )
  

# export histogram

ggsave(
  filename = paste0("Morisita_histograms_", collapse_strategy, "_crab.pdf"),
  plot = g,
  width = 8,
  height = 6
)

ggsave(
  filename = paste0("Morisita_histograms_", collapse_strategy, "_crab.png"),
  plot = g,
  width = 8,
  height = 6
)

saveRDS(g, paste0("Morisita_histograms_", collapse_strategy, "_crab.rds"),)
