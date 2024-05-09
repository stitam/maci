rm(list = ls())

rooted_trees <- readRDS("input/global_ST2_tree/rooted_trees.rds")
tree <- rooted_trees$rtt_rms_1

assemblies <- readRDS("results_redacted/filter_assemblies/aci_filtered.rds")

index <- which(!assemblies$assembly %in% tree$tip.label)
assemblies <- assemblies[-index,]

index <- which(!tree$tip.label %in% assemblies$assembly)
testthat::expect_equal(length(index), 0)

# calculate snps for each rooted tree
snp <- ape::node.depth.edgelength(tree)[1:ape::Ntip(tree)]

# rescale tip_dates to calendar dates
tip_dates <- unname(sapply(tree$tip.label, function(x) {
  assemblies$collection_day[which(assemblies$assembly == x)]
}))

tip_dates <- as.Date(tip_dates, origin = "1970-01-01")

tip_dates <- lubridate::decimal_date(tip_dates)

# recalculate root-to-tip regression using calendar dates
fit <- lm(snp~tip_dates)

# calculate metrics for each fit
results <- data.frame(
  r.squared = summary(fit)$r.squared,
  adj.r.squared = summary(fit)$adj.r.squared,
  rse = summary(fit)$sigma,
  ssr = sum((summary(fit)$residuals)^2),
  mrca = -fit$coef[1]/fit$coef[2],
  first = min(tip_dates, na.rm = TRUE)[1]
)
results$first <- as.Date(results$first, origin = "1970-01-01")
results$mrca <- as.Date(results$mrca, origin = "1970-01-01")

df <- data.frame(
  snp = fit$model$snp,
  date = as.Date(fit$model$tip_dates, origin = "1970-01-01")
)

g <- ggplot(df, aes(date, snp)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Date") +
  ylab("Genetic distance from root (SNP/genome)") +
  theme(
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
  )
 
ggsave(filename = "manual_results/root_to_tip/FigS4A_root_to_tip.pdf")
ggsave(filename = "manual_results/root_to_tip/FigS4A_root_to_tip.png")
