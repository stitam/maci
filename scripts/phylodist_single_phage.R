library(dplyr)
library(ggplot2)
library(ggh4x)
library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-a", "--assemblies"),
    type = "character",
    help = "Path to assemblies file."
  ),
  make_option(
    c("-t", "--trees"),
    type = "character",
    help = "All dated trees collected in a single rds file."
  ),
  make_option(
    c("-s", "--sensitivity"),
    type = "character",
    help = "A table of phage-host sensitivity test results."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    assemblies = "results/filter_assemblies/aci_filtered.rds",
    trees = "data/all_dated_trees.rds",
    sensitivity = "labdata/aci_project/Ab_all_strains_phages_spotassay_PFU.tsv"
  )
}

# import assemblies
assemblies <- readRDS(args$assemblies)

# import trees
trees <- readRDS(args$trees)

# calculate cophenetic (patristic) distance matrix for each pair of genomes
patrdist <- lapply(trees, ape::cophenetic.phylo)
names(patrdist) <- names(trees)

# import sensitivity test results
phres <- read.csv(args$sensitivity, sep = "\t", na.strings = c("", "NA"))

# only keep strains where sequences are available
phres <- phres[which(!is.na(phres$Assembly)),]

# only keep strains where tree for the MLST is available
phres <- phres[which(phres$MLST %in% names(trees)),]

# COMMENT BECAUSE WE WANT TO TEST ALL PHAGES SEPARATELY
# only keep strains which were tested against all phages
# phres <- phres[complete.cases(phres[,8:ncol(phres)]),]

# only keep strains which are present on a tree
# only keep strains which are present on a single tree
remove <- vector()
for (i in phres$Assembly) {
  hit <- FALSE
  for (j in names(trees)) {
    if (i %in% trees[[j]]$tip.label) {
      if (hit == FALSE) {
        hit <- TRUE
      } else {
        msg <- paste0(i, " found on multiple trees!")
        message(msg)
      }
    }
  }
  if (hit == FALSE) {
    msg <- paste0(i, " not found!")
    message(msg)
    remove <- c(remove, i)
  }
}

if (length(remove) > 0) {
  phres <- phres[-which(phres$Assembly %in% remove),]
}

# format phres
for (i in 1:nrow(phres)) {
  for (j in 8:ncol(phres)) {
    # Count partial sensitivity to a phage as no sensitivity
    if (!is.na(phres[i,j]) && phres[i,j] == "1") {
      phres[i,j] <- "0"
    }
    if (!is.na(phres[i,j]) && phres[i,j] %in% c("YES", "Yes")) {
      phres[i,j] <- "1"
    }
    if (!is.na(phres[i,j]) && phres[i,j] == "NO") {
      phres[i,j] <- "0"
    }
  }
}

for (j in 8:ncol(phres)) {
  phres[,j] <- as.numeric(phres[,j])
  phres[,j] <- ifelse(
    phres[,j] > 1,
    1,
    phres[,j]
  )
}

# MODIFY APPLY FUNCTION TO ALLOW NA VALUES
# only keep strains which were sensitive to at least one phage
effective_phage_count <- apply(
  phres[,8:ncol(phres)], 1, function(x) sum(x, na.rm = TRUE))
phres <- phres[which(effective_phage_count > 0),]

# only keep strains that belong to KL types against which at least one phage
# was found. e.g. if there was an effective phage against KL3 then keep all KL3
# strains, regardless of their phage profile.

#infection_count <- sapply(unique(phres$KL), function(x) {
#  sum(phres[which(phres$KL == x),8:ncol(phres)])
#})
#kl_to_keep <- names(infection_count)[which(infection_count > 0)]
#phres <- phres[which(phres$KL %in% kl_to_keep),]

comp <- combn(phres$Assembly, 2) %>% t() %>% as.data.frame()
names(comp) <- c("A1", "A2")
comp$A1_mlst <- sapply(comp$A1, function(x) {
  index <- which(phres$Assembly == x)
  phres$MLST[index]
})
comp$A2_mlst <- sapply(comp$A2, function(x) {
  index <- which(phres$Assembly == x)
  phres$MLST[index]
})
comp$A1_KL <- sapply(comp$A1, function(x) {
  index <- which(phres$Assembly == x)
  phres$KL[index]
})
comp$A2_KL <- sapply(comp$A2, function(x) {
  index <- which(phres$Assembly == x)
  phres$KL[index]
})
comp$A1_region23 <- sapply(comp$A1, function(x) {
  index <- which(assemblies$assembly == x)
  assemblies$region23[index]
})
# Combine Eastern and Southern Europe into a single region
comp$A1_region23 <- ifelse(
  comp$A1_region23 %in% c("eastern_europe", "southern_europe"),
  "east_south_europe",
  comp$A1_region23
)

comp$A2_region23 <- sapply(comp$A2, function(x) {
  index <- which(assemblies$assembly == x)
  assemblies$region23[index]
})

# Combine Eastern and Southern Europe into a single region
comp$A2_region23 <- ifelse(
  comp$A2_region23 %in% c("eastern_europe", "southern_europe"),
  "east_south_europe",
  comp$A2_region23
)

comp$Region <-factor(
  ifelse(
    comp$A1_region23 == comp$A2_region23,
    "Same region",
    "Different regions"
  ), 
  levels = c("Different regions", "Same region")
)
# only compare strains which had the same MLST and K type
comp <- comp[which(comp$A1_mlst == comp$A2_mlst & comp$A1_KL == comp$A2_KL),]

# copy patristic distances from distance matrix
comp$patristic <- NA
for (i in 1:nrow(comp)) {
  ST <- comp$A1_mlst[i]
  index_x <- which(rownames(patrdist[[ST]]) == comp$A1[i])
  index_y <- which(colnames(patrdist[[ST]]) == comp$A2[i])
  comp$patristic[i] <- patrdist[[ST]][index_x, index_y]
}

# since distance is in years, divide by two to relate to divergence times.
comp$patristic <- comp$patristic/2

phagenames <- names(phres)[8:ncol(phres)]

for (j in phagenames) {
  comp[[ncol(comp) + 1]] <- NA
  names(comp)[ncol(comp)] <- j
  for (i in 1:nrow(comp)) {
    index_x <- which(phres$Assembly == comp$A1[i])
    index_y <- which(phres$Assembly == comp$A2[i])
    comp[i, ncol(comp)] <- sum(c(phres[index_x, j], phres[index_y, j]))
  }
}

comp <- tidyr::pivot_longer(
  comp,
  cols = phagenames,
  names_to = "phage",
  values_to = "icount"
)

# define phage-serotype categories
comp$phagesero <- paste0(comp$A1_mlst, "-", comp$A1_KL, " ", comp$phage)

# filter to serotype-phage pairs with at least 10 comparisons (5 strains)
# filter to serotype-phage pairs where at least one strain was sensitive
phagesero <- comp %>% 
  group_by(phagesero) %>% 
  summarise(
    comparison_count = n(),
    isolate_count = (1+sqrt(1+8*comparison_count))/2,
    sensitivity_count = sum(icount, na.rm = TRUE)
  ) %>% 
  filter(
    isolate_count >= 5 &
    sensitivity_count > 0
  )
comp <- comp[which(comp$phagesero %in% phagesero$phagesero),]

# select a representative phage for each scenario and only include those pairs

# convert phage sensitivities to long format
phres_long <- phres %>%
  select(1, 8:ncol(phres)) %>%
  tidyr::pivot_longer(
    cols = 2:ncol(.),
    names_to = "phage",
    values_to = "susceptible"
  )

# merge susceptibilities into comp table
comp$A1_sus = NA
comp$A2_sus = NA
for (i in 1:nrow(comp)) {
  comp$A1_sus[i] <- phres_long$susceptible[which(
    phres_long$Assembly == comp$A1[i] & phres_long$phage == comp$phage[i])]
  comp$A2_sus[i] <- phres_long$susceptible[which(
    phres_long$Assembly == comp$A2[i] & phres_long$phage == comp$phage[i])]
}

# only keep complete rows (otherwise phylodist ranges would be distorted)
comp <- comp[complete.cases(comp),]

reps <- data.frame()
for (i in unique(comp$phagesero)) {
  comp_small <- comp[which(comp$phagesero == i),]
  df <- data.frame(
    phagesero = i,
    sero = strsplit(i, split = " ")[[1]][1],
    phage = strsplit(i, split = " ")[[1]][2],
    rep = unique(c(comp_small$A1, comp_small$A2))
  )
  df$sus = NA
  df$phylodist_from = NA
  df$phylodist_to = NA
  df$phylodist_span = NA
  df$phylodist_sd = NA
  for (j in 1:nrow(df)) {
    df$sus[j] <- phres_long$susceptible[which(
      phres_long$Assembly == df$rep[j] & phres_long$phage == df$phage[j])]
    
    phylorange <- range(comp_small$patristic[which(
      comp_small$A1 == df$rep[j] | comp_small$A2 == df$rep[j]
    )], na.rm = TRUE)
    
    df$phylodist_from[j] <- phylorange[1]
    df$phylodist_to[j] <- phylorange[2]
    df$phylodist_span[j] <- phylorange[2]-phylorange[1]
    
    df$phylodist_sd[j] <- sd(comp_small$patristic[which(
      comp_small$A1 == df$rep[j] | comp_small$A2 == df$rep[j]
    )], na.rm = TRUE)
  }
  reps <- dplyr::bind_rows(
    reps,
    df
  )
}

# set random seed because multiple reps may satisfy filter criteria
set.seed(0)

# reps:
# - are infected by the phage
# - the phyodist range they span with their pairs is min 90% of the max range
# - phylodist sd is maximal
# these are successive filters from top to bottom
reps_final <- reps %>%
  group_by(phagesero) %>%
  filter(sus == 1) %>%
  filter(phylodist_span >= 0.9*max(phylodist_span, na.rm = TRUE)) %>%
  filter(phylodist_sd == max(phylodist_sd)) %>%
  slice_sample(n = 1)

# test that each serotype phage combination has a single rep
testthat::expect_true(nrow(reps_final) == length(unique(reps_final$phagesero)))

comp <- comp %>%
  group_by(phagesero) %>%
  filter(
    A1 == reps_final$rep[which(reps_final$phagesero == unique(phagesero))] | 
    A2 == reps_final$rep[which(reps_final$phagesero == unique(phagesero))]
  )

# prepare plot for each phage-serotype pair
for (i in unique(comp$phagesero)) {
  comp_small <- comp |> filter(phagesero == i)
  if (sum(!is.na(comp_small$icount)) < 10) next()
  if (sum(comp_small$icount, na.rm = TRUE) == 0) next()
  g <- ggplot(comp_small, aes(patristic, icount)) + 
    geom_point(aes(col = Region)) + 
    theme_bw() +
    xlab("Average time of divergence (years)") +
    ylab("Number of isolates successfully infected from a pair of isolates") +
    ggtitle(paste0(
      i, " Number of isolates: ",
      length(unique(c(comp_small$A1, comp_small$A2)))
    ))
  
  #used for testing
  #ggsave(
  #  filename = paste0(i, ".pdf"),
  #  plot = g
  #)  
}

# filter to selected phage-serotype pairs
comp <- comp[which(comp$phagesero %in% c(
  "ST2-KL2 Rocket",
  "ST2-KL2 Fishpie",
  "ST2-KL2 Silvergun",
  "ST2-KL3 Rocket",
  "ST2-KL3 Highwayman",
  "ST2-KL3 Silvergun",
  "ST636-KL40 Rocket",
  "ST636-KL40 PhT2.v2",
  "ST636-KL40 Konradin.v2"
)),]

# rename phagesero to include desired strip text
for (i in unique(comp$phagesero)) {
  index <- which(comp$phagesero == i)
  count <- length(unique(c(
    comp$A1[which(comp$phagesero == i)],
    comp$A2[which(comp$phagesero == i)]
  )))
  comp$phagesero[index] <- paste0(
    comp$phagesero[index], "\n",
    "N = ",
    count
  )
}

new_phagesero <- unique(comp$phagesero)

# order phagesero factor levels for plotting
comp$phagesero <- factor(comp$phagesero, levels = c(
  new_phagesero[grep("ST2-KL2 Rocket", new_phagesero)],
  new_phagesero[grep("ST2-KL2 Fishpie", new_phagesero)],
  new_phagesero[grep("ST2-KL2 Silvergun", new_phagesero)],
  new_phagesero[grep("ST2-KL3 Rocket", new_phagesero)],
  new_phagesero[grep("ST2-KL3 Highwayman", new_phagesero)],
  new_phagesero[grep("ST2-KL3 Silvergun", new_phagesero)],
  new_phagesero[grep("ST636-KL40 Rocket", new_phagesero)],
  new_phagesero[grep("ST636-KL40 PhT2.v2", new_phagesero)],
  new_phagesero[grep("ST636-KL40 Konradin.v2", new_phagesero)]
))

# prepare composite plot for selected phage-serotype pairs

strip <- ggh4x::strip_themed(
  background_x = elem_list_rect(fill = c(
    rep("grey", times = 1),
    rep("#00A36C", times = 2),
    rep("grey", times = 1),
    rep("#00A36C", times = 2),
    rep("grey", times = 3)
  ))
)

#Use 0-1 scale instead of 1-2 scale
comp$icount <- comp$icount - 1

g <- ggplot(comp, aes(patristic, icount)) + 
  geom_point(
    position = position_jitter(height = 0.05, seed = 0),
    alpha = 0.5,
    aes(col = Region)
  ) +
  geom_smooth(method = "loess", span = 0.85, size = 0.1) +
  facet_wrap2(~phagesero, ncol = 3, strip = strip) +
  theme_bw() +
  xlab("Average time of divergence\nsince common ancestor with sensitive reference isolate (years)") +
  ylab("Probability of successful phage infection") +
  theme(
    plot.title = ggtext::element_markdown(size = 8),
    plot.margin = unit(c(0,0,0,0), "cm"),
    panel.background = element_rect(fill = "white"),
    axis.line=element_line(color="black", linewidth=0.1),
    axis.ticks = element_line(linewidth=0.1),
    axis.text = element_text(size = 8),
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom"
  )

g

ggsave(
  filename = "phage_serotype_sensitivity.pdf",
  plot = g,
  width = 15,
  height = 15,
  units = "cm"
)

ggsave(
  filename = "phage_serotype_sensitivity.png",
  plot = g,
  width = 15,
  height = 15,
  units = "cm"
)
