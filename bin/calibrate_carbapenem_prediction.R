library(devtools)
library(dplyr)
library(ggplot2)
library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to the project directory."
  ),
  make_option(
    c("-f", "--file"),
    type = "character",
    help = "Path to aci prediction results."
  ),
  make_option(
    c("-g", "--resgene_groups"),
    type = "character",
    help = "Path to a reference table containing '-like' classifications"
  ),
  make_option(
    c("-a", "--amr_db"),
    type = "character",
    help = "Path to a tidy amr data set."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "aci",
    file = "results/run_qc_checks/aci_with_qc.rds",
    resgene_groups = "aci/data/resgene_groups.csv",
    amr_db = "results/tidy_bvbrc/amr_bvbrc.rds"
  )
}

load_all(args$project_dir)
aci <- readRDS(args$file)
resgene_groups <- read.csv(args$resgene_groups, sep = "\t")
amr_db <- readRDS(args$amr_db)

# prepare input variables

aci <- aci[which(aci$qc_pass == TRUE), ]

# define new variable that contains resgene groups instead of resgenes
aci$beta_lactam_groups <- sapply(aci$beta_lactam, function(x) {
  genes <- strsplit(x, "\\|")[[1]]
  for (i in seq_along(genes)) {
    index <- which(resgene_groups$gene == genes[i])
    if (length(index) == 1) {
      genes[i] <- resgene_groups$group[index]
    }
  }
  out <- paste(genes, collapse = "|")
  return(out)
})

# convert aci table to long format using the beta-lactam groups
aci_wide <- aci[,which(names(aci) %in% c(
  "biosample",
  "ISAba1",
  "carO",
  "beta_lactam_groups"
))]
aci_wide <- aci_wide[complete.cases(aci_wide),]
aci_wide <- tidyr::separate_longer_delim(
  data = aci_wide,
  cols = beta_lactam_groups,
  delim = "|"
)

if (any(aci_wide$beta_lactam_groups == "NA")) {
  aci_wide <- aci_wide[-which(aci_wide$beta_lactam_groups == "NA"),]
}

aci_wide <- tidyr::pivot_wider(
  aci_wide,
  id_cols = c("biosample", "ISAba1", "carO"),
  names_from = beta_lactam_groups,
  values_from = beta_lactam_groups,
  values_fn = unique
)

# do not touch the first three columns
for (j in 4:ncol(aci_wide)) {
  aci_wide[,j] <- !is.na(aci_wide[,j])
}

saveRDS(aci_wide, "aci_wide.rds")

# prepare output variables

carbapenems <- amr_db$Antibiotic[grep("penem$", amr_db$Antibiotic)] %>%
  sort() %>%
  unique()

# the column "carbapenem" seems to be a subset of "imipenem". Remove.
carbapenems <- carbapenems[-which(carbapenems == "carbapenem")]

amr_db <- amr_db[which(amr_db$Antibiotic %in% carbapenems),]

amr_wide <- tidyr::pivot_wider(
  amr_db,
  id_cols = "BioSample.Accession",
  names_from = "Antibiotic",
  values_from = "Resistant.Phenotype",
  values_fn = function(x) paste(sort(unique(x)), collapse = "|")
)

# define function to force binary classification
# TODO add reference that intermediate > resistant.
rename_amr <- function(x) {
  x <- ifelse(x == "Intermediate", "Resistant", x)
  x <- ifelse(x == "Intermediate|Resistant", "Resistant", x)
  x <- ifelse(x == "Intermediate|Susceptible", "Susceptible", x)
  x <- ifelse(x == "Intermediate|Unknown", "Resistant", x)
  x <- ifelse(x == "Resistant|Unknown", "Resistant", x)
  x <- ifelse(x == "Susceptible|Unknown", "Susceptible", x)
  x <- ifelse(x == "Resistant|Susceptible", NA, x)
  x <- ifelse(x == "Intermediate|Resistant|Susceptible", NA, x)
  return(x)
}

# force binary classification for carbapenems
for (j in 2:ncol(amr_wide)) {
  amr_wide[,j] <- sapply(amr_wide[,j], rename_amr)
}

# check that resistance is "Resistant", "Susceptible" or NA.
test <- amr_wide[,2:ncol(amr_wide)] %>% unlist() %>% sort() %>% unique()
testthat::expect_true(all(test %in% c("Resistant", "Susceptible")))

# define new variable whether isolate is resistant against at least one
# carbapenem. NOTE: if an isolate is resistant against one antibiotic and
# susceptible against another, this will return TRUE > This is a permissive
# classifier.
amr_wide$resmin1 <- sapply(1:nrow(amr_wide), function(x) {
    profile <- vector()
    for (i in carbapenems) {
        profile <- c(profile, amr_wide[[which(names(amr_wide) == i)]][x])
    }
    if ("Resistant" %in% profile) return("Resistant")
    if ("Susceptible" %in% profile) return("Susceptible")
    return(NA)
})

# define new variable whether isolate is resistant against all carbapenems
# where a measurement result is available. NOTE: if an isolate is resistant
# against one antibiotic and susceptible against another, this will return 
# FALSE > This is a restrictive classifier.
amr_wide$res_all <- sapply(1:nrow(amr_wide), function(x) {
    profile <- vector()
    for (i in carbapenems) {
        profile <- c(profile, amr_wide[[which(names(amr_wide) == i)]][x])
    }
    if ("Susceptible" %in% profile) return("Susceptible")
    if ("Resistant" %in% profile) return("Resistant")
    return(NA)
})

saveRDS(amr_wide, "amr_wide.rds")

# Analysis with resmin1

df <- dplyr::right_join(
  aci_wide,
  amr_wide[,which(names(amr_wide) %in% c("BioSample.Accession", "resmin1"))], 
  by = c("biosample" = "BioSample.Accession")
)

# remove entries which do not have assemblies 
df <- df[complete.cases(df),]

index_all_true <- which(apply(df[2:ncol(df)], 2, function(x) {
  all(x == TRUE, na.rm = TRUE)
}))

if (length(index_all_true) > 0) {
  index_adjusted <- index_all_true + 1
  df <- df[-index_adjusted]
}

index_all_false <- which(apply(df[2:ncol(df)], 2, function(x) {
  all(x == FALSE, na.rm = TRUE)
}))

if (length(index_all_false) > 0) {
  index_adjusted <- index_all_false + 1
  df <- df[-index_adjusted]
}

# tbl: table(real, predicted)
classification_stats <- function(tbl) {
  precision <- tbl[1,2]/(tbl[1,2] + tbl[2,2])
  tpr <- tbl[1,2]/(tbl[1,2] + tbl[1,1])
  fpr <- tbl[2,2]/(tbl[2,2] + tbl[2,1])
  return(list(
    "table" = tbl,
    "Precision" = signif(precision, 3),
    "Recall" = signif(tpr, 3),
    "FPR" = signif(fpr, 3)
  ))
}

# stats for individual features

dets <- matrix(NA, nrow = length(2:(ncol(df)-1)), ncol = 9)
colnames(dets) <- c("determinant", "FN", "TN", "TP", "FP", "p_value", "Precision", "Recall", "FPR")
dets[,1] <- names(df)[2:(ncol(df)-1)]

# save table for contingency table calculation
df_contingency <- df

for (i in 1:nrow(dets)) {
  tbl <- table(df$resmin1, df[[dets[i,1]]])
  f <- try(fisher.test(df$resmin1, df[[dets[i,1]]]), silent = TRUE)
  if (!inherits(f, "try-error")){
    fp <- signif(f$p.value, 3)
  } else {
    fp = NA
  }
  cl <- classification_stats(tbl)
  dets[i, 2:ncol(dets)] <- c(
    cl$table[1,1],
    cl$table[2,1],
    cl$table[1,2],
    cl$table[2,2],
    fp,
    cl$Precision,
    cl$Recall,
    cl$FPR
  )
}

dets <- as.data.frame(dets)

# order entries by decreasing precision
# super important for multiple determinants analysis!
dets <- dets[order(dets$Precision, decreasing = TRUE),]

dets$p_value_fdr <- p.adjust(dets$p_value, method = "fdr")
dets <- dplyr::relocate(dets, p_value_fdr, .after = "p_value")

saveRDS(dets, file = "determinants.rds")

write.table(
  dets,
  file = "determinants.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# for multiple determinants analysis, filter dets table to significant entries
dets <- dets[which(dets$p_value_fdr < 0.05),]

multidets <- matrix(NA, nrow = length(2:(ncol(df)-1)), ncol = 9)
colnames(multidets) <- c("count", "FN", "TN", "TP", "FP", "p_value", "Precision", "Recall", "FPR")

dets_old <- dets

# Remove ISAba1 because it is not potential resistance gene by itself
dets <- dets[-which(dets$determinant == "ISAba1"),]

# start with predictor with highest precision
# keep adding the predictor with the next highest precision
for (k in 1:nrow(dets)) {
  multidets[k,1] <- k
  predictors <- dets[1:k,1]
  predicted <- rep(FALSE, times = nrow(df))
  for (i in 1:nrow(df)) {
    if (any(df[i, which(names(df) %in% predictors)] == TRUE)) {
      predicted[i] <- TRUE
    }
  }
  if (all(predicted) == TRUE) break()
  tbl <- table(df$resmin1, predicted)
  f <- try(fisher.test(df$resmin1, predicted), silent = TRUE)
  if (!inherits(f, "try-error")){
    fp <- signif(f$p.value, 3)
  } else {
    fp = NA
  }
  cl <- classification_stats(tbl)
  multidets[k, 2:ncol(multidets)] <- c(
    cl$table[1,1],
    cl$table[2,1],
    cl$table[1,2],
    cl$table[2,2],
    fp,
    cl$Precision,
    cl$Recall,
    cl$FPR
  )
}

multidets <- as.data.frame(multidets)
multidets <- multidets[complete.cases(multidets),]

saveRDS(multidets, file = "multiple_determinants.rds")

precision_recall <- ggplot(multidets, aes(Recall, Precision)) + 
  geom_point() +
  xlim(c(0,1)) +
  ylim(c(0,1))

ggsave(
  filename = "precision_recall.pdf",
  plot = precision_recall,
  width = 4,
  height = 4
)

roc <- ggplot(multidets, aes(FPR, Recall)) + 
  geom_point() +
  xlim(c(0,1)) +
  ylim(c(0,1))

ggsave(
  filename = "roc.pdf",
  plot = roc,
  width = 4,
  height = 4
)

# look at ISAba1

isaba <- matrix(NA, nrow = length(2:(ncol(df)-1)), ncol = 9)
colnames(isaba) <- c("determinant", "FN", "TN", "TP", "FP", "p_value", "Precision", "Recall", "FPR")
isaba[,1] <- names(df)[2:(ncol(df)-1)]

for (k in 1:nrow(isaba)) {
  predictors <- c("ISAba1", unname(isaba[k,1]))
  predicted <- rep(FALSE, times = nrow(df))
  for (i in 1:nrow(df)) {
    if (all(df[i, which(names(df) %in% predictors)] == TRUE)) {
      predicted[i] <- TRUE
    }
  }
  if (all(predicted == TRUE) == TRUE) next()
  if (all(predicted == FALSE) == TRUE) next()
  tbl <- table(df$resmin1, predicted)
  f <- try(fisher.test(df$resmin1, predicted), silent = TRUE)
  if (!inherits(f, "try-error")){
    fp <- signif(f$p.value, 3)
  } else {
    fp = NA
  }
  cl <- classification_stats(tbl)
  isaba[k, 2:ncol(isaba)] <- c(
    cl$table[1,1],
    cl$table[2,1],
    cl$table[1,2],
    cl$table[2,2],
    fp,
    cl$Precision,
    cl$Recall,
    cl$FPR
  )
}

saveRDS(isaba, file = "ISAba1.rds")

# prep contingency table for the carbapenem definition that was selected manually.

df2 <- df_contingency

# based on names(df2) these are the genes that are present at least once and
# are also listed in predict_carbapenem.R
df2$predicted <- df2$`blaOXA-23-like` |
  df2$`blaOXA-40-like` |
  df2$`blaOXA-58-like` |
  df2$`blaOXA-143-like` |
  df2$`blaOXA-312` |
  df2$`blaGES-11` |
  df2$`blaIMP-16` |
  df2$`blaNDM-1`

tbl <- table(df2$resmin1, df2$predicted)
f <- try(fisher.test(df2$resmin1, df2$predicted), silent = TRUE)
if (!inherits(f, "try-error")){
  fp <- signif(f$p.value, 3)
} else {
  fp = NA
}
cl <- classification_stats(tbl)

