library(devtools)
library(ggplot2)
library(ggtext)
library(optparse)
library(patchwork)
rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to the project directory."
  ),
  make_option(
    c("-a", "--aci_path"),
    type = "character",
    help = "Path to aci prediction results."
  ),
  make_option(
    c("-m", "--minyear_threshold"),
    type = "integer",
    help = "Lowest year to include."
  ),
  make_option(
    c("-c", "--mincount_threshold"),
    type = "integer",
    help = "Lowest count to include."
  ),
  make_option(
    c("-s", "--selected_serotypes"),
    type = "integer",
    help = "Path to a table containing a list of selected serotypes"
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "aci",
    aci_path = "results/filter_assemblies/aci_filtered.rds",
    minyear_threshold = 2009,
    mincount_threshold = 50,
    selected_serotypes = "results/calc_serotype_freqs/geodate/global_or_prevalent_serotypes_region23_collapse_geodate.tsv"
  )
}

load_all(args$project_dir)

minyear_threshold <- as.numeric(args$minyear_threshold)

# import data set
aci <- readRDS(args$aci_path)

# import list of selected serotypes
prevalent <- read.csv(args$selected_serotypes, sep = "\t")

aci <- aci[, c(
  "assembly",
  "mlst",
  "k_serotype",
  "serotype",
  "continent",
  "region23",
  "country",
  "collection_year"
)]

# make data set wider
nc <- ncol(aci)
aci <- tidyr::pivot_wider(
  aci,
  names_from = serotype,
  values_from = serotype
)
aci[, nc: ncol(aci)] <- apply(
  aci[, nc: ncol(aci)], 2, function(x) !is.na(x)
)

#initialise logistic regression summary tibble
logres <- tibble::tibble()

# Look at each serotype
# Does its frequency increase worldwide, corrected for continents?
# Does its frequency increase continent wide, corrected for regions?
# Count for each serotype the number of continents where year is significant

for (i in nc: ncol(aci)) {
  index <- which(aci[[i]] == TRUE & aci$collection_year >= minyear_threshold)
  if (length(index) > 1) {
    # 'minyear' is the year of the first isolation since 'minyear_threshold' 
    minyear <-min(aci$collection_year[
      which(aci[[i]] == TRUE & aci$collection_year >= minyear_threshold)])
    # logistic regression worldwide
    # continent is added to predictors
    aci_crop <- aci[which(aci$collection_year >= minyear), ]
    if (length(unique(aci_crop$collection_year)) > 1) {
      if (length(unique(aci_crop$continent)) > 1) {
        fit_log <- suppressWarnings(glm(
          aci_crop[[names(aci_crop)[i]]] ~  aci_crop$collection_year + aci_crop$continent,
          family = "binomial"
        ))
        new_res <- tibble::tibble(
          serotype = names(aci_crop)[i],
          scope = "world",
          count = sum(aci_crop[[names(aci_crop)[i]]]),
          by = "collection_year, continent",
          predictor = row.names(summary(fit_log)$coefficients),
          estimate = summary(fit_log)$coefficients[,1],
          p_value = summary(fit_log)$coefficients[,4]
        )
        new_res <- new_res[-which(new_res$predictor == "(Intercept)"),]
        new_res$predictor <- gsub("aci_crop\\$", "", new_res$predictor)
        new_res$predictor <- gsub("continent", "", new_res$predictor)
      } else {
        fit_log <- suppressWarnings(glm(
          aci_crop[[names(aci_crop)[i]]] ~  aci_crop$collection_year,
          family = "binomial"
        ))
        new_res <- tibble::tibble(
          serotype = names(aci_crop)[i],
          scope = "world",
          count = sum(aci_crop[[names(aci_crop)[i]]]),
          by = "collection_year",
          predictor = row.names(summary(fit_log)$coefficients),
          estimate = summary(fit_log)$coefficients[,1],
          p_value = summary(fit_log)$coefficients[,4]
        )
        new_res <- new_res[-which(new_res$predictor == "(Intercept)"),]
        new_res$predictor <- gsub("aci_crop\\$", "", new_res$predictor)
        new_res$predictor <- gsub("continent", "", new_res$predictor)
      }
    logres <- dplyr::bind_rows(logres, new_res)  
    }
    for (j in unique(aci$continent)) {
      index2 <- which(aci[[i]] == TRUE & 
                       aci$continent == j & 
                       aci$collection_year >= minyear_threshold)
      # only progress if a serotype was found at least once
      if (length(index2) > 1) {
        minyear2 <-min(aci$collection_year[index2])
        aci_crop2 <- aci[which(aci$continent == j & aci$collection_year >= minyear2), ]
        if (length(unique(aci_crop2$collection_year)) > 1) {
          if (length(unique(aci_crop2$region23)) > 1) {
            fit_log <- suppressWarnings(glm(
              aci_crop2[[names(aci_crop2)[i]]] ~  aci_crop2$collection_year + aci_crop2$region23,
              family = "binomial"
          ))
          new_res <- tibble::tibble(
            serotype = names(aci_crop)[i],
            scope = j,
            count = sum(aci_crop2[[names(aci_crop2)[i]]]),
            by = "collection_year, region23",
            predictor = row.names(summary(fit_log)$coefficients),
            estimate = summary(fit_log)$coefficients[,1],
            p_value = summary(fit_log)$coefficients[,4]
          )
          new_res <- new_res[-which(new_res$predictor == "(Intercept)"),]
          new_res$predictor <- gsub("aci_crop2\\$", "", new_res$predictor)
          new_res$predictor <- gsub("region23", "", new_res$predictor)
          } else {
            fit_log <- suppressWarnings(glm(
              aci_crop2[[names(aci_crop2)[i]]] ~  aci_crop2$collection_year,
              family = "binomial"
          ))
          new_res <- tibble::tibble(
            serotype = names(aci_crop)[i],
            scope = j,
            count = sum(aci_crop2[[names(aci_crop2)[i]]]),
            by = "collection_year",
            predictor = row.names(summary(fit_log)$coefficients),
            estimate = summary(fit_log)$coefficients[,1],
            p_value = summary(fit_log)$coefficients[,4]
          )
          new_res <- new_res[-which(new_res$predictor == "(Intercept)"),]
          new_res$predictor <- gsub("aci_crop2\\$", "", new_res$predictor)
          new_res$predictor <- gsub("region23", "", new_res$predictor)
          }
          logres <- dplyr::bind_rows(logres, new_res)
        }
      }
    }
  }
}

logres$estimate <- signif(logres$estimate, 2)
logres$p_value <- signif(logres$p_value, 5)

# export logistic regression results
saveRDS(logres, file = "logistic_regression_results.rds")
write.table(
  logres,
  file = "logistic_regression_results.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# adjust p values for collection year, due to multiple comparison
logres_date <- logres[which(
  logres$predictor == "collection_year" & 
  logres$serotype %in% prevalent$serotype &
  logres$count >= args$mincount_threshold),]
logres_date$p_value_corr <- stats::p.adjust(logres_date$p_value, method = "fdr")
logres_date$p_value_corr <- signif(logres_date$p_value_corr, 3)

# export results
saveRDS(logres_date, file = "logistic_regression_results_with_fdr.rds")
write.table(
  logres_date,
  file = "TableS3_logistic_regression_results_with_fdr.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# tidy logistic regression results
index <- which(logres_date$p_value_corr < 0.05)
logres2 <- logres_date[index,c(
  "serotype", "scope", "estimate"
)]
logres2$estimate <- signif(logres2$estimate, 3)
logres3 <- tidyr::pivot_wider(
  logres2,
  names_from = scope,
  values_from = estimate
)
logres3 <- dplyr::right_join(
  dplyr::distinct(logres[which(logres$scope == "world"), c("serotype","count")]),
  logres3,
  by = "serotype"
)
logres3$nup <- apply(logres3[, 4:ncol(logres3)], 1, function(x) {
  length(which(sign(x) == 1))
})
logres3$ndown <- apply(logres3[, 4:ncol(logres3)], 1, function(x) {
  length(which(sign(x) == -1))
})

# export tidy results
saveRDS(logres3, file = "logistic_regression_results_tidy.rds")
write.table(
  logres3,
  file = "logistic_regression_results_tidy.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
)

# create figures
figdir <- "figures"
if (!dir.exists(figdir)) dir.create(figdir)
# ratio of isolates over time - barplot
composite_fig <- list()

for (i in 1:nrow(logres_date)) {
  sero <- logres_date$serotype[i]
  if (logres_date$scope[i] == "world") {
    aci_crop <- aci
    index <- which(
      aci[[sero]] == TRUE & 
      aci$collection_year >= minyear_threshold
    )
  } else {
    aci_crop <- aci[which(aci$continent == logres_date$scope[i]),]
    index <- which(
      aci[[sero]] == TRUE &
      aci$continent == logres_date$scope[i] &
      aci$collection_year >= minyear_threshold
    )
  }
  minyear <-min(aci$collection_year[index])
  #ggsave(
  #  filename = paste0(
  #    figdir, "/ratio_of_",
  #    sero, "_isolates_by_year_", logres_date$scope[i], ".png"),
  #  plot = plot_ts(
  #    aci_crop,
  #    var = sero,
  #    type = "ratio",
  #    geom = "geom_col",
  #    show_count = TRUE,
  #    window = 1,
  #    minyear = minyear
  #  ),
  #  width = 8,
  #  height= 6
  #) 
  if (logres_date$p_value_corr[i] <= 0.05 & logres_date$scope[i] != "world") {
    figname <- paste0(
      gsub("_", "-", sero), " ", 
      Hmisc::capitalize(logres_date$scope[i])
    )
    figtitle <- paste0(
      figname,
      "<br>",
      "p = ", signif(logres_date$p_value_corr[i], 3),
      rep(" ", times = 1)
    )
    if (sign(logres_date$estimate[i]) == 1) {
      figtitle <- paste0(
        figtitle,
        "<span style = 'color: red;font-weight:bold;'>\u2b06\u2b06\u2b06</span>"
      )
    }
    if (sign(logres_date$estimate[i]) == -1) {
      figtitle <- paste0(
        figtitle,
        "<span style = 'color: forestgreen;font-weight:bold;'>\u2b07\u2b07\u2b07</span>"
      )
    }
    ts_df <- get_ts_df(
      df = aci_crop,
      var = sero,
      timevar = "collection_year", 
      window = 1
    )
    ts_df <- ts_df[which(
      ts_df$collection_year >= minyear &
      ts_df[[sero]] == TRUE
    ),]
    ts_df$collection_year <- factor(
      as.character(ts_df$collection_year),
      levels = sort(unique(ts_df$collection_year))
    )
    composite_fig[[figname]] <- ggplot(ts_df, aes(collection_year, ratio_window)) +
      geom_col(fill = "#A9A9A9", alpha = 0.9) +
      labs(
        title = figtitle,
        x = "",
        y = "") +
      annotate(
        "text",
        label = round(ts_df$yearcount, 0),
        x = ts_df$collection_year,
        y = max(ts_df$ratio_window),
        size = 8/.pt,
        hjust = 1,
        angle = 90
      ) + 
      theme(
        plot.title = ggtext::element_markdown(size = 8),
        plot.margin = unit(c(0,0,0,0), "cm"),
        panel.background = element_rect(fill = "white"),
        axis.line=element_line(color="black", linewidth=0.1),
        axis.ticks = element_line(linewidth=0.1),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90),
      )
  }
}

nc <- ifelse(length(composite_fig) >=3, 3, length(composite_fig))
nr <- ceiling(length(composite_fig) / nc)

cfig <- ggpubr::ggarrange(
  plotlist = composite_fig,
  legend = "none",
  ncol = nc,
  nrow = nr
)
cfig <- ggpubr::annotate_figure(
  cfig,
  left = grid::textGrob(
    "Relative prevalence",
    rot = 90,
    vjust = 1,
    gp = grid::gpar(fontsize = 8)
  ),
  bottom = grid::textGrob(
    "Year of isolation",
    gp = grid::gpar(fontsize = 8)
  )
)

ggsave(
  filename = "figures/FigS23_ratio_of_serotypes_by_year_composite.png",
  plot = cfig,
  width = 18,
  units = "cm"
)

ggsave(
  filename = "figures/FigS23_ratio_of_serotypes_by_year_composite.pdf",
  plot = cfig,
  width = 18,
  units = "cm",
  device = cairo_pdf
)



# # ratio of isolates over time - lines
# for (i in 1:nrow(logres_date)) {
#   sero <- logres_date$serotype[i]
#   minyear <-min(aci$collection_year[
#     which(aci[[sero]] == TRUE & aci$collection_year >= minyear_threshold)])
#   if (logres_date$scope[i] == "world") {
#     aci_crop <- aci
#     ggsave(
#       filename = paste0(
#         figdir, "/ratio_of_",
#         sero, "_isolates_by_year_", logres_date$scope[i], "_lines.pdf"),
#       plot = plot_ts(
#         aci_crop,
#         var = sero,
#         by = "continent",
#         type = "ratio",
#         window = 2,
#         minyear = minyear),
#       width = 10,
#       height= 8
#     )
#   } else {
#     aci_crop <- aci[which(aci$continent == logres_date$scope[i]),]
#     ggsave(
#       filename = paste0(
#         figdir, "/ratio_of_",
#         sero, "_isolates_by_year_", logres_date$scope[i], "_lines.pdf"),
#       plot = plot_ts(
#         aci_crop,
#         var = sero,
#         by = "region23",
#         type = "ratio",
#         window = 2,
#         minyear = minyear),
#       width = 10,
#       height= 8
#     )
#   }
# }
