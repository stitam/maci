library(dplyr)
library(ggplot2)
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
    assemblies = "results_redacted/filter_assemblies/aci_filtered.rds",
    trees = "data/all_dated_trees.rds",
    sensitivity = "data/phage_sensitivity_measurements/TableS2_T10.tsv"
  )
}

# import assemblies
assemblies <- readRDS(args$assemblies) %>% dplyr::filter(filtered)

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

# only keep strains which were tested against all phages
phres <- phres[complete.cases(phres[,7:ncol(phres)]),]

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
  for (j in 7:ncol(phres)) {
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

for (j in 7:ncol(phres)) {
  phres[,j] <- as.numeric(phres[,j])
  phres[,j] <- ifelse(
    phres[,j] > 1,
    1,
    phres[,j]
  )
}

# only keep strains which were sensitive to at least one phage
effective_phage_count <- apply(phres[,7:ncol(phres)], 1, sum)
phres <- phres[which(effective_phage_count > 0),]

# only keep strains that belong to KL types against which at least one phage
# was found. e.g. if there was an effective phage against KL3 then keep all KL3
# strains, regardless of their phage profile.

#infection_count <- sapply(unique(phres$KL), function(x) {
#  sum(phres[which(phres$KL == x),7:ncol(phres)])
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
  comp$A1_region23 %in% c("Eastern Europe", "Southern Europe"),
  "east_south_europe",
  comp$A1_region23
)

comp$A2_region23 <- sapply(comp$A2, function(x) {
  index <- which(assemblies$assembly == x)
  assemblies$region23[index]
})

# Combine Eastern and Southern Europe into a single region
comp$A2_region23 <- ifelse(
  comp$A2_region23 %in% c("Eastern Europe", "Southern Europe"),
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

# calculate phage sensitivity distance matrix
phagedist <- as.matrix(vegan::vegdist(phres[7:ncol(phres)], "jaccard", na.rm = FALSE))
rownames(phagedist) <- phres$Assembly
colnames(phagedist) <- phres$Assembly

# copy phage sensitivity distances from distance matrix
comp$phagedist <- NA
for (i in 1:nrow(comp)) {
  index_x <- which(rownames(phagedist) == comp$A1[i])
  index_y <- which(colnames(phagedist) == comp$A2[i])
  comp$phagedist[i] <- phagedist[index_x, index_y]
}

comp$patristic <- signif(comp$patristic, 4)
comp$phagedist <- signif(comp$phagedist, 4)

# Some phage distances may be NaN depending on the distance metric. Remove these.
comp <- comp[complete.cases(comp),]

write.table(
  comp,
  file = "phylodist_phagedist.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# define combined serotype
phres$serotype <- paste0(phres$MLST, " - ", phres$KL)
comp$serotype <- paste0(comp$A1_mlst, " - ", comp$A1_KL)

# filter to serotypes with at least 45 comparisons (10 strains)
serotypes <- comp %>% 
  group_by(serotype) %>% 
  summarise(count = n()) %>% 
  filter(count >= 45)
comp <- comp[which(comp$serotype %in% serotypes$serotype),]

# PERFORM MANTEL TEST
mantel <- list()
for (i in unique(comp$serotype)) {
  ST <- strsplit(i, " - ")[[1]][1]
  assembly_sub <- phres$Assembly[which(phres$serotype == i)]
  # PATRISTIC DISTANCES
  # select the distance matrix for the ST
  patrdist_sub <- patrdist[[ST]]
  # subset the distance matrix to ST-K serotype
  index_x <- which(rownames(patrdist_sub) %in% assembly_sub)
  index_y <- which(colnames(patrdist_sub) %in% assembly_sub)
  patrdist_sub <- patrdist_sub[index_x, index_y]
  # sort by row and column names
  patrdist_sub <- patrdist_sub[order(rownames(patrdist_sub)),]
  patrdist_sub <- patrdist_sub[,order(colnames(patrdist_sub))]
  # PHAGE DISTANCES
  index_x <- which(rownames(phagedist) %in% assembly_sub)
  index_y <- which(colnames(phagedist) %in% assembly_sub)
  phagedist_sub <- phagedist[index_x, index_y]
  # sort by row and column names
  phagedist_sub <- phagedist_sub[order(rownames(phagedist_sub)),]
  phagedist_sub <- phagedist_sub[,order(colnames(phagedist_sub))]
  # CONSISTENCY CHECKS
  testthat::expect_true(all(rownames(patrdist_sub) == colnames(patrdist_sub)))
  testthat::expect_true(all(rownames(phagedist_sub) == colnames(phagedist_sub)))
  testthat::expect_true(all(rownames(patrdist_sub) == rownames(phagedist_sub)))
  testthat::expect_true(all(colnames(patrdist_sub) == colnames(phagedist_sub)))
  # convert to dist
  patrdist_sub <- as.dist(patrdist_sub, diag = TRUE)
  phagedist_sub <- as.dist(phagedist_sub, diag = TRUE)
  # perform mantel test
  if (any(is.nan(phagedist_sub))) next()
  mantel[[i]] <- vegan::mantel(patrdist_sub, phagedist_sub, method = "spearman")
}
comp$mantel_spearman <- sapply(comp$serotype, function(x) {
  if (x %in% names(mantel)) {
    mantel[[x]]$signif
  } else {
    return(NA)
  }
})

# Prepare Plot

comp_jitter <- comp

set.seed(0)
comp_jitter$patristic <- jitter(comp$patristic, amount = 1)
comp_jitter$phagedist <- jitter(comp$phagedist, amount = 0.03)

# define limits for plotting
xmin <- min(comp_jitter$patristic, na.rm = TRUE)
xmax <- max(comp_jitter$patristic, na.rm = TRUE)
ymin <- min(comp_jitter$phagedist, na.rm = TRUE)
ymax <- max(comp_jitter$phagedist, na.rm = TRUE)

g <- ggplot(comp_jitter, aes(patristic, phagedist)) + 
  geom_point(alpha = 0.5, aes(col = Region), size  = 0.2) +
  geom_smooth(method = "loess", span = 5, size = 0.1) +
  xlim(xmin, xmax) +
  ylim(ymin, ymax) +
  xlab("Average time of divergence (years)") + 
  ylab("Dissimilarity of phage susceptibility profiles \n (Jaccard distance)") +
  facet_wrap(serotype~., nrow = 1)

# add mantel test p values
mantel_df <- data.frame(
  patristic = 22,
  phagedist = 0.1,
  serotype = unique(comp$serotype)
)
mantel_df$mantel_spearman <- sapply(mantel_df$serotype, function(x) {
  paste0(
    "Mantel test significance: \n",
    if (x %in% names(mantel)) {
      mantel[[x]]$signif
    } else {
      return(NA)
    }
  )})

# g <- g + geom_text(data = mantel_df, aes(label = mantel_spearman), size = 2)

ggsave(
  filename = "phylodist_phagedist.pdf",
  plot = g,
  units = "cm",
  height = 15,
  width = 30
)

ggsave(
  filename = "phylodist_phagedist.png",
  plot = g,
  units = "cm",
  height = 15,
  width = 30
)

saveRDS(g, "phylodist_phagedist.rds")

# ANALYSE PAIRS WHICH HAD IDENTICAL PHAGE SENSITIVITY PROFILES

# Assign strains with identical phage sensitivity profiles to the same cluster
d <- dist(phres[,7:ncol(phres)])
h <- hclust(d)
phres$cluster <- cutree(h, h = 0)

# tidy cluster names
nchar_max <- max(nchar(phres$cluster))
phres$cluster <- sapply(phres$cluster, function(x) {
  x <- as.character(x)
  n <- nchar(x)
  zeros <- nchar_max - n
  if(zeros == 0) {
    out2 <- paste0("cluster_", x)
  }
  if (zeros > 0) {
    out <- strsplit(x, "")[[1]]
    out2 <- rep(NA, times = nchar_max + 8)
    out2[1:8] <- c("c","l","u","s","t","e","r","_")
    out2[9:(9+zeros-1)] <- "0"
    out2[(9+zeros):length(out2)] <- out
    out2 <- paste(out2, collapse = "")
  }
  return(out2)
})

write.table(
  phres,
  file = "phage_host_resistance_profiles.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
