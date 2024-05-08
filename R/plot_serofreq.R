#' Plot frequency of serotypes in data set
#' 
#' This function visualises the diversity of Multi Locus Sequence Types (MLSTs),
#' K Locus types and predicted resistance genes within the studied data set. 
#' @param df a data fram returned by the nextflow pipeline
#' @param min_count integer; the lowest ST-KL serogroup count to include
#' @param max_bin_count integer, the maximum number of resistance family bins
#' to include
#' @import ggplot2
#' @importFrom reshape2 dcast
#' @importFrom tidyr pivot_longer
#' @export
plot_serofreq <- function(df,
                          min_count = 1,
                          max_bin_count = 3) {
  df$mlst_count <- sapply(df$mlst, function(x) length(which(df$mlst == x)))
  df$k_serotype_count <- sapply(df$k_serotype,
                                function(x) length(which(df$k_serotype == x)))
  df$serogroup_count <- mapply(function(x, y) {
    length(which(df$k_serotype==x & df$mlst == y))
  }, df$k_serotype, df$mlst)
  resclasses <- c(
    "aminoglycoside",
    "beta_lactam",
    "colistin",
    "disinfectant",
    "fosfomycin",
    "fusidicacid",
    "glycopeptide",
    "macrolide",
    "nitroimidazole",
    "oxazolidinone",
    "phenicol",
    "pseudomonicacid",
    "quinolone",
    "rifampicin",
    "sulphonamide",
    "tetracycline",
    "trimethoprim"
  )
  df$res_family_count <- apply(
    df[,which(names(df) %in% resclasses)], 1, function(x) sum(!is.na(x)))
  df$serogroup_res_family_count_median <- mapply(function(x,y) {
    index <- which(df$k_serotype==x & df$mlst == y)
    median(df$res_family_count[index], na.rm = TRUE)
  }, df$k_serotype, df$mlst)
  c <- max(df$serogroup_res_family_count_median, na.rm = TRUE)
  get_breaks <- function(c, nbreaks){
    round(seq(from = 0, to = ceiling(c), by = 1/nbreaks*c),0)
  }
  df$serogroup_res_family_count_median <- cut(
    df$serogroup_res_family_count_median,
    breaks = get_breaks(c, max_bin_count))
  df_dcast <- reshape2::dcast(
    df, k_serotype ~ mlst, value.var = "assembly", fun.aggregate = length)
  df_longer <- tidyr::pivot_longer(
    df_dcast, !k_serotype, names_to = "mlst", values_to = "count")
  index <- which(df$serogroup_count >= min_count)
  mlst_to_plot <- unique(df$mlst[index])
  k_serotype_to_plot <- unique(df$k_serotype[index])
  index <- which(
    df$mlst %in% mlst_to_plot & df$k_serotype %in% k_serotype_to_plot)
  ggplot(df[index, ], aes(mlst, k_serotype))+
    geom_point(aes(
      size = serogroup_count,
      col = serogroup_res_family_count_median))+
    theme(axis.text.x = element_text(angle = 90))+
    xlab("MLST") +
    ylab("K serotype")+
    guides(
      size = guide_legend(title = "serogroup count"),
      col = guide_legend(title = "resistance families (median)"))
}
