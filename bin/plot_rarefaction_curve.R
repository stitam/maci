library(dplyr)
library(ggplot2)
library(optparse)
library(scales)
library(vegan)

rm(list = ls())

# create log file and start logging
con <- file(paste0("log.txt"))
sink(con, split = TRUE)

args_list <- list(
  make_option(
    c("-f", "--file"),
    type = "character",
    help = "Path to aci prediction results."
  ),
  make_option(
    c("-m", "--minyear"),
    type = "character",
    help = "Minimum year to plot."
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
    file = "aci_collapse_geodate_crab.rds",
    minyear = "2016",
    regions = "geographic_locations_in_study.tsv"
  )
}

collapse_strategy <- args$file %>% strsplit(split = "_") %>% unlist()
collapse_strategy <- collapse_strategy[3]

# import data
aci <- readRDS(args$file)

# filter by collection_year
aci <- aci[which(aci$collection_year >= as.numeric(args$minyear)),]

# import region pretty names and colors
regions <- read.csv(args$regions, sep = "\t", na = "") %>%
  dplyr::filter(variable == "region23") %>%
  dplyr::select(term, term_pretty, colour) %>%
  dplyr::rename(region23 = term, region23_pretty = term_pretty)

region_colors <- regions$colour
names(region_colors) <- regions$region23_pretty

# make region names pretty for plotting
aci$region23 <- unname(sapply(aci$region23, function(x) {
  index <- which(regions$region23 == x)
  if (length(index) == 0) {
    msg <- paste0("Region not found: ", x)
    stop(msg)
    return(x)
  } else {
    return(regions$region23_pretty[index])
  }
}))

# import country pretty names and colors
countries <- read.csv(args$regions, sep = "\t", na = "") %>%
  dplyr::filter(variable == "country") %>%
  dplyr::select(term, term_pretty, colour) %>%
  dplyr::rename(country = term, country_pretty = term_pretty)

country_colors <- countries$colour
names(country_colors) <- countries$country_pretty

# make country names pretty for plotting
aci$country <- unname(sapply(aci$country, function(x) {
  index <- which(countries$country == x)
  if (length(index) == 0) {
    msg <- paste0("Country not found: ", x)
    stop(msg)
    return(x)
  } else {
    return(countries$country_pretty[index])
  }
}))

# analyse regions

aci_wide <- aci %>% 
  select(assembly, region23, serotype) %>%
  tidyr::pivot_wider(
    id_cols = region23,
    names_from = serotype,
    values_from = serotype,
    values_fn = function(x) length(!is.na(x))
  )

# replace NAs with 0s
for (i in 1:nrow(aci_wide)) {
  for (j in 1:ncol(aci_wide)) {
    if (is.na(aci_wide[i, j])) {
      aci_wide[i, j] <- 0
    }
  }
}

aci_wide_rn <- aci_wide %>% tibble::column_to_rownames("region23")

plot_df <- vegan::rarecurve(aci_wide_rn, tidy = TRUE)

LT_names <- as.character(unique(plot_df$Site))
LT <- ifelse(LT_names %in% c("Northern America", "Eastern Asia"), "dashed", "solid")
names(LT) <- LT_names

# consistency check
# check that all regions have colors
testthat::expect_true(
  all(plot_df$Site %in% names(region_colors)))

# main plot
xmax <- max(plot_df$Sample[which(!plot_df$Site %in% c(
  "Northern America", "Eastern Asia"
))])
ymax <- max(plot_df$Species[which(!plot_df$Site %in% c(
  "Northern America", "Eastern Asia"
))])
g1 <- ggplot(plot_df, aes(Sample, Species)) + 
  geom_line(aes(col = Site, linetype = Site)) + 
  scale_colour_manual(values = region_colors, name="Region") +
  scale_linetype_manual(values = LT, name = "Region")+
  xlim(0, xmax) +
  ylim(0, ymax) +
  ylab("Number of MLST-CPS types")

# embedded plot
g2 <- ggplot(plot_df, aes(Sample, Species)) + 
  geom_line(aes(col = Site, linetype = Site)) + 
  scale_colour_manual(values = region_colors, name="Region") +
  scale_linetype_manual(values = LT, name = "Region")+
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = scales::alpha(colour = "white", alpha = 0.5),),
    panel.background = element_rect(fill = "transparent"),
    panel.grid = element_blank()
  ) +
  xlab("") +
  ylab("")

g3 <- g1 + 
  annotation_custom(
    ggplotGrob(g2),
    xmin = 0.5 * xmax,
    xmax = xmax,
    ymin = 0,
    ymax = 0.5 * ymax
  )

print(sort(table(aci$region23), descending = TRUE))

ggsave(
  filename = paste0(
    "FigS2A_rarecurve_crab_collapse_",
    collapse_strategy,
    "_minyear_",
    args$minyear, "_regions.pdf"
  ),
  plot = g3,
  width = 10,
  height = 6
)
ggsave(
  filename = paste0(
    "FigS2A_rarecurve_crab_collapse_",
    collapse_strategy,
    "_minyear_",
    args$minyear, "_regions.png"
  ),
  plot = g3,
  width = 10,
  height = 6
)
saveRDS(
  object = g3,
  file = paste0(
    "FigS2A_rarecurve_crab_collapse_",
    collapse_strategy,
    "_minyear_",
    args$minyear, "_regions.rds"
  )
)

# analyse European countries

aci_wide <- aci %>% 
  filter(continent == "europe") %>%
  select(assembly, country, serotype) %>%
  tidyr::pivot_wider(
    id_cols = country,
    names_from = serotype,
    values_from = serotype,
    values_fn = function(x) length(!is.na(x))
  )

# replace NAs with 0s
for (i in 1:nrow(aci_wide)) {
  for (j in 1:ncol(aci_wide)) {
    if (is.na(aci_wide[i, j])) {
      aci_wide[i, j] <- 0
    }
  }
}

aci_wide_rn <- aci_wide %>% tibble::column_to_rownames("country")

plot_df <- vegan::rarecurve(aci_wide_rn, tidy = TRUE)

LT_names <- as.character(unique(plot_df$Site))
LT <- ifelse(LT_names %in% c("China", "USA"), "dashed", "solid")
names(LT) <- LT_names

# consistency check
# check that all regions have colors
testthat::expect_true(
  all(plot_df$Site %in% names(country_colors)))

# main plot
# xmax <- max(plot_df$Sample[which(!plot_df$Site %in% c(
#   "China", "USA"
# ))])
# ymax <- max(plot_df$Species[which(!plot_df$Site %in% c(
#   "China", "USA"
# ))])
g1 <- ggplot(plot_df, aes(Sample, Species)) + 
  geom_line(aes(col = Site, linetype = Site)) + 
  scale_colour_manual(values = country_colors, name="Country") +
  scale_linetype_manual(values = LT, name = "Country") +
  guides(col = guide_legend(title="Country")) +
  # xlim(0, xmax) +
  # ylim(0, ymax) +
  ylab("Number of MLST-CPS types")

# embedded plot
# g2 <- ggplot(plot_df, aes(Sample, Species)) + 
#   geom_line(aes(col = Site, linetype = Site)) + 
#   scale_colour_manual(values = country_colors, name="Country") +
#   scale_linetype_manual(values = LT, name = "Country")+
#   theme(
#     legend.position = "none",
#     plot.background = element_rect(fill = scales::alpha(colour = "white", alpha = 0.5),),
#     panel.background = element_rect(fill = "transparent"),
#     panel.grid = element_blank()
#   )+
#   xlab("") +
#   ylab("")

# g3 <- g1 + 
#   annotation_custom(
#     ggplotGrob(g2),
#     xmin = 0.5 * xmax,
#     xmax = xmax,
#     ymin = 0,
#     ymax = 0.5 * ymax
#   )

# g3

print(sort(table(aci$country), descending = TRUE))

ggsave(
  filename = paste0(
    "FigS2B_rarecurve_Europe_crab_collapse_",
    collapse_strategy,
    "_minyear_",
    args$minyear, "_countries.pdf"
  ),
  plot = g1,
  width = 10,
  height = 6
)
ggsave(
  filename = paste0(
    "FigS2B_rarecurve_Europe_crab_collapse_",
    collapse_strategy,
    "_minyear_",
    args$minyear, "_countries.png"
  ),
  plot = g1,
  width = 10,
  height = 6
)

# analyse non-European countries

aci_wide <- aci %>% 
  filter(continent != "europe") %>%
  select(assembly, country, serotype) %>%
  tidyr::pivot_wider(
    id_cols = country,
    names_from = serotype,
    values_from = serotype,
    values_fn = function(x) length(!is.na(x))
  )

# replace NAs with 0s
for (i in 1:nrow(aci_wide)) {
  for (j in 1:ncol(aci_wide)) {
    if (is.na(aci_wide[i, j])) {
      aci_wide[i, j] <- 0
    }
  }
}

aci_wide_rn <- aci_wide %>% tibble::column_to_rownames("country")

plot_df <- vegan::rarecurve(aci_wide_rn, tidy = TRUE)

LT_names <- as.character(unique(plot_df$Site))
LT <- ifelse(LT_names %in% c("China", "USA"), "dashed", "solid")
names(LT) <- LT_names

# consistency check
# check that all regions have colors
testthat::expect_true(
  all(plot_df$Site %in% names(country_colors)))

# main plot
xmax <- max(plot_df$Sample[which(!plot_df$Site %in% c(
  "China", "USA"
))])
ymax <- max(plot_df$Species[which(!plot_df$Site %in% c(
  "China", "USA"
))])
g1 <- ggplot(plot_df, aes(Sample, Species)) + 
  geom_line(aes(col = Site, linetype = Site)) + 
  scale_colour_manual(values = country_colors, name="Country") +
  scale_linetype_manual(values = LT, name = "Country") +
  guides(col = guide_legend(title="Country")) +
  xlim(0, xmax) +
  ylim(0, ymax) +
  ylab("Number of MLST-CPS types")

# embedded plot
g2 <- ggplot(plot_df, aes(Sample, Species)) + 
  geom_line(aes(col = Site, linetype = Site)) + 
  scale_colour_manual(values = country_colors, name="Country") +
  scale_linetype_manual(values = LT, name = "Country")+
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = scales::alpha(colour = "white", alpha = 0.5),),
    panel.background = element_rect(fill = "transparent"),
    panel.grid = element_blank()
  )+
  xlab("") +
  ylab("")

g3 <- g1 + 
  annotation_custom(
    ggplotGrob(g2),
    xmin = 0.5 * xmax,
    xmax = xmax,
    ymin = 0,
    ymax = 0.5 * ymax
  )

g3

print(sort(table(aci$country), descending = TRUE))

ggsave(
  filename = paste0(
    "FigS2B2_rarecurve_not_Europe_crab_collapse_",
    collapse_strategy,
    "_minyear_",
    args$minyear, "_countries.pdf"
  ),
  plot = g3,
  width = 10,
  height = 6
)
ggsave(
  filename = paste0(
    "FigS2B2_rarecurve_not_Europe_crab_collapse_",
    collapse_strategy,
    "_minyear_",
    args$minyear, "_countries.png"
  ),
  plot = g3,
  width = 10,
  height = 6
)

# # analyse states in usa

# # only keep entries where the country is the USA and the state is not NA
# aci_usa <- aci[which(aci$country == "USA" & !is.na(aci$gadm_level1)), ]

# aci_wide <- aci_usa %>% 
#   select(assembly, gadm_level1, serotype) %>%
#   tidyr::pivot_wider(
#     id_cols = gadm_level1,
#     names_from = serotype,
#     values_from = serotype,
#     values_fn = function(x) length(!is.na(x))
#   )

# # replace NAs with 0s
# for (i in 1:nrow(aci_wide)) {
#   for (j in 1:ncol(aci_wide)) {
#     if (is.na(aci_wide[i, j])) {
#       aci_wide[i, j] <- 0
#     }
#   }
# }

# aci_wide_rn <- aci_wide %>% tibble::column_to_rownames("gadm_level1")

# plot_df <- vegan::rarecurve(aci_wide_rn, tidy = TRUE)

# g4 <- ggplot(plot_df, aes(Sample, Species)) + 
#   geom_line(aes(col = Site)) + 
#   guides(col = guide_legend(title="State"))
#   ylab("Number of ST-K serotypes")

# print(sort(table(aci_usa$gadm_level1), descending = TRUE))

# ggsave(
#   filename = paste0(
#     "FigSX_rarecurve_crab_collapse_",
#     collapse_strategy,
#     "_minyear_",
#     args$minyear, "_usa.pdf"
#   ),
#   plot = g4,
#   width = 10,
#   height = 6
# )
# ggsave(
#   filename = paste0(
#     "FigSX_rarecurve_crab_collapse_",
#     collapse_strategy,
#     "_minyear_",
#     args$minyear, "_usa.png"
#   ),
#   plot = g4,
#   width = 10,
#   height = 6
# )
# saveRDS(
#   object = g4,
#   file = paste0(
#     "FigSX_rarecurve_crab_collapse_",
#     collapse_strategy,
#     "_minyear_",
#     args$minyear, "_usa.rds"
#   )
# )

# # analyse provinces in china

# # only keep entries where the country is China and the province is not NA
# aci_china <- aci[which(aci$country == "China" & !is.na(aci$gadm_level1)), ]

# aci_wide <- aci_china %>% 
#   select(assembly, gadm_level1, serotype) %>%
#   tidyr::pivot_wider(
#     id_cols = gadm_level1,
#     names_from = serotype,
#     values_from = serotype,
#     values_fn = function(x) length(!is.na(x))
#   )

# # replace NAs with 0s
# for (i in 1:nrow(aci_wide)) {
#   for (j in 1:ncol(aci_wide)) {
#     if (is.na(aci_wide[i, j])) {
#       aci_wide[i, j] <- 0
#     }
#   }
# }

# aci_wide_rn <- aci_wide %>% tibble::column_to_rownames("gadm_level1")

# plot_df <- vegan::rarecurve(aci_wide_rn, tidy = TRUE)

# g5 <- ggplot(plot_df, aes(Sample, Species)) + 
#   geom_line(aes(col = Site)) + 
#   guides(col = guide_legend(title="Province"))
#   ylab("Number of ST-K serotypes")

# print(sort(table(aci_china$gadm_level1), descending = TRUE))

# ggsave(
#   filename = paste0(
#     "FigSX_rarecurve_crab_collapse_",
#     collapse_strategy,
#     "_minyear_",
#     args$minyear, "_china.pdf"
#   ),
#   plot = g5,
#   width = 10,
#   height = 6
# )

# ggsave(
#   filename = paste0(
#     "FigSX_rarecurve_crab_collapse_",
#     collapse_strategy,
#     "_minyear_",
#     args$minyear, "_china.png"
#   ),
#   plot = g5,
#   width = 10,
#   height = 6
# )

# saveRDS(
#   object = g5,
#   file = paste0(
#     "FigSX_rarecurve_crab_collapse_",
#     collapse_strategy,
#     "_minyear_",
#     args$minyear, "_china.rds"
#   )
# )

# end logging
sink(con)
