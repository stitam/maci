library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    "--tree",
    type = "character",
    help = "Path to the tree file"
  ),
  make_option(
    "--tree_tbl",
    type = "character",
    help = "Path to the tree table file"
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    tree = "data/global_ST2_tree/dated_tree.nwk",
    tree_tbl = "data/global_ST2_tree/tree_tbl.rds"
  )
}

library(dplyr)
library(ggplot2)

# import files
tree <- ape::read.tree(args$tree)
tree_tbl <- readRDS(args$tree_tbl)

# define focal countries
# TODO this should not be hard-coded
focal_countries <- c(
  "Bosnia and Herzegovina",
  "Hungary",
  "Montenegro",
  "Romania",
  "Serbia"
)

testthat::expect_true(all(focal_countries %in% tree_tbl$country))

# tips that belong to focal countries
focal_tips <- tree_tbl$node[which(tree_tbl$country %in% focal_countries)]
focal_tips <- focal_tips[which(focal_tips <= length(tree$tip.label))]

# capital names
capitals <- c(
  "belgrade",
  "bucharest",
  "city11",
  "podgorica"
)

expand_geo <- function(x, capitals) {
  if (grepl(" ", x)) {
    out <- strsplit(x, " ")[[1]]
    out = list(
      city = ifelse(
        out[1] %in% capitals,
        "capital",
        "not_capital"
      ),
      country = out[2],
      original_category = x
    )
    out$new_category <- paste(out$city, out$country, sep = " ")
  } else {
    out <-  list (
      city = NA,
      country = x,
      original_category = x,
      new_category = x
    )
  }
  return(out)
}

# assemble data frame, one row for each tip

df <- data.frame(
  focal_tip = focal_tips,
  tip_date = sapply(focal_tips, function(x) {
    tree_tbl$collection_day[which(tree_tbl$node == x)]
  }),
  tip_city = sapply(focal_tips, function(x) {
    tree_tbl$city_pooled[which(tree_tbl$node == x)]
  })
)

df$tip_date <- as.Date(df$tip_date, origin = "1970-01-01")

df$tip_city_aggregated <- sapply(df$tip_city, function(x) {
  expand_geo(x, capitals = capitals)$new_category
})

df$tip_country <- sapply(df$focal_tip, function(x) {
  tree_tbl$country[which(tree_tbl$node == x)]
})

# This is the node "BEFORE!!" which the state change may occur
# This node has the same character state as the query node
get_closest_sc_node <- function(tree_tbl, node) {
  sc_prob <- 0
  node_i <- node
  while(sc_prob == 0) {
    sc_prob <- tree_tbl$city_pooled_sc_prob[which(tree_tbl$node == node_i)]
    if (sc_prob == 0) {
      node_i <- tree_tbl$parent[which(tree_tbl$node == node_i)]
    }
  }
  return(node_i)
}

df$closest_sc_node <- sapply(df$focal_tip, function(x) {
  get_closest_sc_node(tree_tbl = tree_tbl, node = x)
})

df$tips_from_node <- sapply(df$closest_sc_node, function(x) {
  df %>% dplyr::filter(closest_sc_node == x) %>% nrow()
})

df$singeton <- df$tips_from_node == 1

df$sc_node_date <- sapply(df$closest_sc_node, function(x) {
  tree_tbl$collection_day[which(tree_tbl$node == x)]
})
df$sc_node_date <- as.Date(df$sc_node_date, origin = "1970-01-01")

df$chain_length <- round(as.numeric(df$tip_date-df$sc_node_date)/365.25,3)

df$from_geo <- sapply(df$closest_sc_node, function(x) {
  index <- which(tree_tbl$node == x)
  tree_tbl$city_pooled_from[index]
})

df$domestic <- NA
df$between_focal <- NA
df$non_focal <- NA

for (i in 1:nrow(df)) {
  from_vec <- strsplit(df$from_geo[i], "\\|")[[1]]
  from_list <- lapply(from_vec, function(x) {
    expand_geo(x, capitals = capitals)
  })
  to <- expand_geo(df$tip_city[i], capitals = capitals)
  intro_type <- vector()
  for (j in 1:length(from_list)) {
    if (is.na(to$city)) {
      next()
    } else if (from_list[[j]]$country == to$country) {
      intro_type <- c(intro_type, "domestic")
    } else if (from_list[[j]]$country %in% c(
      "(BA)", "(HU)", "(ME)", "(RO)", "(RS)"
    )) {
      intro_type <- c(intro_type, "between_focal")
    } else if (from_list[[j]]$original_category %in% c(
      "other_eastern_european_country",
      "other_southern_european_country",
      "other_european_country",
      "other_country_outside_europe")) {
      intro_type <- c(intro_type, "non_focal")
    } else {
      stop()
    }
  }
  df$domestic[i] <- round(length(which(intro_type == "domestic"))/length(from_vec), 3)
  df$between_focal[i] <- round(length(which(intro_type == "between_focal"))/length(from_vec), 3)
  df$non_focal[i] <- round(length(which(intro_type == "non_focal"))/length(from_vec), 3)
}

# assemble data frame, one row for each transmission lineage

dfsum <- df %>%
  filter(tips_from_node > 1) %>%
  group_by(closest_sc_node) %>%
  summarise(
    city = unique(tip_city),
    city_aggregated = unique(tip_city_aggregated),
    country = unique(tip_country),
    longest_chain = max(chain_length),
    from_geo = unique(from_geo),
    domestic = max(domestic),
    between_focal = max(between_focal),
    non_focal = max(non_focal),
    .groups = "drop"
  )

dfsum$longest_chain <- NA
for (i in seq_along(dfsum$closest_sc_node)) {
  index <- which(df$closest_sc_node == dfsum$closest_sc_node[i])
  dfsum$longest_chain[i] <- max(df$chain_length[index])
}

# assign each transmission to either "Domestic", "Between focal" or "Non-focal"

dfsum$group <- NA
set.seed(0)
for (i in 1:nrow(dfsum)) {
  dfsum$group[i] <- sample(
    c("Domestic", "Between focal", "Non-focal"),
    size = 1,
    prob = c(dfsum$domestic[i], dfsum$between_focal[i], dfsum$non_focal[i])
  )
}

dfsum$group <- factor(
  dfsum$group, levels =  c("Domestic", "Between focal", "Non-focal")
)

# remove unnecessary columns
dfsum <- dfsum %>%
  select(-city_aggregated, -domestic, -between_focal, -non_focal)

# make column names more informative

dfsum <- dfsum %>%
  dplyr::rename(
    `MRCA node` = closest_sc_node,
    `location of lineage (city)` = city,
    `location of lineage (country)` = country,
    `geographical source` = from_geo,
    `duration of lineage (years)` = longest_chain,
    `type of lineage` = group
  )

write.table(
  df,
  file = "chain_length_from_closest_transmission.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  dfsum,
  file = "longest_chain_from_closest_transmission_no_singletons.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

dfsum <- dfsum %>%
  group_by(`type of lineage`) %>%
  mutate(group_mean = mean(`duration of lineage (years)`))

# PLOT 1
# - transmission lineage duration, dotplot

g1 <- ggplot(dfsum, aes(x = `type of lineage`, y = `duration of lineage (years)`)) + 
  geom_dotplot(binaxis = "y", stackdir = "center") +
  geom_errorbar(
    aes(ymin = group_mean, ymax = group_mean),
    linewidth = 0.1,
    width = 0.5,
    col = "red")+
  ylab("Duration of transmission\nlineage (years)") +
  xlab("Source of transmission lineage") +
  guides(fill = "none") +
  theme(
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.y=element_line(color="gray", linetype = 3, linewidth = 0.3),
    axis.line=element_line(color="black", linewidth=0.1),
    axis.ticks = element_line(linewidth=0.1)
  )

ggsave(
  filename = "transmission_lineages.pdf",
  plot = g1,
  width = 6,
  height = 2
)
ggsave(
  filename = "transmission_lineages.png",
  plot = g1,
  width = 6,
  height = 2
)
saveRDS(g1, file = "transmission_lineages.rds")

# # pie chart, should be embedded into plot 1 (?)

# introsum <- dfsum %>%
#   group_by(group) %>%
#   summarise(count = n())
# introsum$percentage <- round(100 * introsum$count/sum(introsum$count), 1)

# pie <- ggplot(introsum, aes(x="", y=percentage, fill=group)) +
#   geom_bar(stat="identity", width=1, color = "white") +
#   coord_polar("y", start=0) +
#   theme_void()
