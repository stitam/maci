#' Check quality of assembled genomes
#' 
#' This function performs a number of deterministic and stochastic QC checks to
#' identify assemblies that may be of poor quality. Deterministic checks compare
#' the value of a variable with a fixed threshold, stochastic checks identify
#' statistically extreme values.
#' @param df data.frame; an assembly summary table.
#' @param taxid character; The NCBI Taxonomy ID for the organism. See Details
#' for more information.
#' @param new_na character; what should the function do if validation introduces
#' new NAs? Either \code{"ignore"} or \code{"stop"}.
#' @param verbose logical; Should verbose messages be printed to the console?
#' @return A list with three elements:
#' \itemize{
#'   \item \code{"data"}: input data frame with new variable \code{qc_pass}. The
#'   new variable is \code{TRUE} if the assembly passed all QC checks, otherwise
#'   \code{FALSE}.
#'   \item \code{"failed"}: data frame of failed checks
#'   \item \code{"thresholds"}: a data frame of stochastic parameter thresholds
#' }
#' @details
#' Some QC checks will run for each organism but some checks will only run for
#' selected organisms. For example, QC checks around capsular types will only
#' run for Acinetobacter baumannii or Klebsiella pneumoniae. This is because 
#' other organisms may not have capsules, or there may not be a typing scheme
#' for the capsule, or it is not yet implemented. Currently the following
#' \code{taxid} values have dedicated tests: \code{470} (Acinetobacter 
#' baumannii), \code{573} (Klebsiella pneumoniae).
#' @details
#' Deterministic checks: \code{"mlst"} is not \code{NA}; exactly 7 genes
#' detected (470); \code{"k_serotype"} is not \code{NA} (470, 573); 
#' \code{"k_confidence"} is at least \code{"Good"} (470, 573); 
#' \code{"kraken2_taxid"} contains the taxid; \code{"coverage"} is at least 25; 
#' BUSCO complete score is at least 95%. Deterministic quality filters are 
#' applied once.
#' @details
#' Stochastic checks: \code{"gc_content"} (use IQR extreme outliers, keep 
#' middle); \code{"contig_count"} (use IQR extreme outliers, keep lower); 
#'  \code{"longest_contig"} (use IQR extreme outliers, keep higher);
#'  \code{"N50"} (use IQR extreme outliers, keep higher); \code{"N95"} (use IQR
#'  extreme outliers, keep higher); \code{"n_count"} (use IQR extreme outliers, 
#'  keep lower, exclude 0); \code{"kraken2_taxid_freq"} (use IQR extreme
#'  outliers, keep higher). Stochastic quality checks are applied iteratively
#'  until no more rows are marked as outliers.
#' @examples
#' \dontrun{
#' check_qc(df, taxid = 470, verbose = TRUE)
#' }
check_qc <- function(
    df,
    taxid, 
    new_na = "ignore",
    verbose = getOption(verbose)
    ) {
  
  if (!"data.frame" %in% class(df)) {
    stop("'df' must be a data frame.")
  }
  
  taxid <- as.numeric(taxid)
  
  if (length(taxid) > 1) {
    stop("'taxid' must be a single integer not a vector or integers.")
  }
  
  df <- validate_typing(df, taxid = taxid, new_na = new_na, verbose = verbose)
  
  if ("qc_pass" %in% names(df)) {
    df <- dplyr::select(df, -qc_pass)
  }

  df_failed <- data.frame()
  
  # add new variable for qc
  df$qc_pass <- TRUE
  
  print("Number of rows:")
  print(nrow(df))
  
  # DETERMINISTIC QUALITY FILTERS
  
  # mlst
  
  print(names(df))
  
  print("Keep rows where MLST is not NA.")
  
  if (length(which(is.na(df$mlst))) > 0) {
    df_rm <- data.frame(
      assembly = df$assembly[which(is.na(df$mlst))],
      reason = "MLST is NA"
    )
    df_failed <- dplyr::bind_rows(
      df_failed,
      df_rm
    )
    df$qc_pass[which(df$assembly %in% df_rm$assembly)] <- FALSE
  }
  
  print("QC passed overall:")
  print(sum(df$qc_pass))
  
  if (taxid == 470) {
    
    print("Keep rows where exactly 7 MLST genes are detected.")
    if (length(which(df$mlst_genes_detected != 7)) > 0) {
      df_rm <- data.frame(
        assembly = df$assembly[which(df$mlst_genes_detected != 7)],
        reason = "More or less than 7 MLST genes detected"
      )
      df_failed <- dplyr::bind_rows(
        df_failed,
        df_rm
      )
      df$qc_pass[which(df$assembly %in% df_rm$assembly)] <- FALSE
    }
    
    print("QC passed overall:")
    print(sum(df$qc_pass))
  }
  
  # k_serotype
  
  if (taxid %in% c(470, 573)) {
    
    print("Keep rows where K serotype is not NA.")
    if (length(which(is.na(df$k_serotype))) > 0) {
      df_rm <- data.frame(
        assembly = df$assembly[which(is.na(df$k_serotype))],
        reason = "K serotype is NA"
      )
      df_failed <- dplyr::bind_rows(
        df_failed,
        df_rm
      )
      df$qc_pass[which(df$assembly %in% df_rm$assembly)] <- FALSE
    }
    
    print("QC passed overall:")
    print(sum(df$qc_pass))
    
    print("Keep rows where K confidence is at least 'Good'.")
    if (length(which(df$k_confidence %in% c("None", "Low"))) > 0) {
      df_rm <- data.frame(
        assembly = df$assembly[which(df$k_confidence %in% c("None", "Low"))],
        reason = "K confidence is less than 'Good'"
      )
      df_failed <- dplyr::bind_rows(
        df_failed,
        df_rm
      )
      df$qc_pass[which(df$assembly %in% df_rm$assembly)] <- FALSE
    }
    
    print("QC passed overall:")
    print(sum(df$qc_pass))
  }
  
  # kraken2_taxid
  
  print(paste0("Keep rows where kraken2_taxid contains '", taxid, "'."))
  index_keep <- strsplit(df$kraken2_taxid, "\\|") %>%
    sapply(function(x) taxid %in% x) %>%
    which()
  
  if (length(which(!index_keep)) > 0) {
    df_rm <- data.frame(
      assembly = df$assembly[which(!index_keep)],
      reason = paste0("kraken2_taxid does not contain '", taxid, "'")
    )
    df_failed <- dplyr::bind_rows(
      df_failed,
      df_rm
    )
    df$qc_pass[which(df$assembly %in% df_rm$assembly)] <- FALSE
  }
  
  print("QC passed overall:")
  print(sum(df$qc_pass))
  
  # coverage
  
  print("Keep rows where coverage is at least 25")
  
  if (length(which(df$coverage < 25)) > 0) {
    df_rm <- data.frame(
      assembly = df$assembly[which(df$coverage < 25)],
      reason = "Coverage is too low"
    )
    df_failed <- dplyr::bind_rows(
      df_failed,
      df_rm
    )
    df$qc_pass[which(df$assembly %in% df_rm$assembly)] <- FALSE
  }
  
  print("QC passed overall:")
  print(sum(df$qc_pass))
  
  # busco
  
  print("Keep rows where BUSCO Complete score is at least 95")
  
  if (length(which(df$busco_complete < 95)) > 0) {
    df_rm <- data.frame(
      assembly = df$assembly[which(df$busco_complete < 95)],
      reason = "BUSCO Complete score is too low"
    )
    df_failed <- dplyr::bind_rows(
      df_failed,
      df_rm
    )
    df$qc_pass[which(df$assembly %in% df_rm$assembly)] <- FALSE
  }
  
  print("QC passed overall:")
  print(sum(df$qc_pass))
  
  # STOCHASTIC QUALITY FILTERS
  
  # IMPORTANT: THE THRESHOLDS BELOW ARE CALCULATED USING ISOLATES THAT HAVE
  # PASSED THE DETERMINISTIC QUALITY FILTERS ABOVE. THIS IS TO ENSURE THAT THE
  # THRESHOLDS ARE NOT BIASED BY ISOLATES THAT ARE CLEARLY OF POOR QUALITY.
  
  ip0 <- 1
  ip1 <- 0
  counter <- 0
  while (ip1 < ip0 & counter < 10) {
    
    print(paste0("Iteration: ", counter))
    
    ip0 <- sum(df$qc_pass)
    
    df_qc <- df[which(df$qc_pass), ]
    
    gc_content <- flag_outliers(
      df_qc$gc_content,
      qlow = 0.25,
      qhigh = 0.75,
      m = 3,
      keep = "middle"
    )
    
    contig_count <- flag_outliers(
      df_qc$contig_count,
      qlow = 0.25,
      qhigh = 0.75,
      m = 3,
      keep = "lower"
    )
    
    if (!"longest_contig" %in% names(df)) {
      stop("'df' must have a variable called 'longest_contig'.")
    }
    
    # TODO: Add assertion that longest contig must be numeric (or integer).
    
    longest_contig <- flag_outliers(
      df_qc$longest_contig,
      qlow = 0.25,
      qhigh = 0.75,
      m = 3,
      keep = "higher"
    )
    
    if (!"N50" %in% names(df)) {
      stop("'df' must have a variable called 'N50'.")
    }
    
    # TODO: Add assertion that N50 must be numeric (or integer).
    
    N50 <- flag_outliers(
      df_qc$N50,
      qlow = 0.25,
      qhigh = 0.75,
      m = 3,
      keep = "higher"
    )
    
    if (!"N95" %in% names(df)) {
      stop("'df' must have a variable called 'N95'.")
    }
    
    # TODO: Add assertion that N95 must be numeric (or integer).
    
    N95 <- flag_outliers(
      df_qc$N95,
      qlow = 0.25,
      qhigh = 0.75,
      m = 3,
      keep = "higher"
    )
    
    n_count <- flag_outliers(
      df_qc$n_count,
      qlow = 0.25,
      qhigh = 0.75,
      m = 3,
      keep = "lower",
      exclude = 0
    )
    
    if (!"kraken2_taxid_freq" %in% names(df)) {
      stop("'df' must have a variable called 'kraken2_taxid_freq'.")
    }
    if (!"numeric" %in% class(df[["kraken2_taxid_freq"]])) {
      stop("'kraken2_taxid_freq' variable within 'df' must be a 'numeric'.")
    }
    
    kraken2_taxid_freq <- flag_outliers(
      df_qc$kraken2_taxid_freq,
      qlow = 0.25,
      qhigh = 0.75,
      m = 3,
      keep = "higher"
    )
    
    qc <- data.frame(
      assembly = df_qc$assembly,
      gc_content = gc_content$keep,
      contig_count = contig_count$keep,
      longest_contig = longest_contig$keep,
      N50 = N50$keep,
      N95 = N95$keep,
      n_count = n_count$keep,
      kraken2_taxid_freq = kraken2_taxid_freq$keep
    )
    
    print(apply(qc[,-1], 2, table))
    
    thresholds <- dplyr::bind_rows(
      gc_content = data.frame(
        name = "gc_content",
        qlow = gc_content$range[1],
        qhigh = gc_content$range[2]
      ),
      contig_count = data.frame(
        name = "contig_count",
        qlow = contig_count$range[1],
        qhigh = contig_count$range[2]
      ),
      longest_contig = data.frame(
        name = "longest_contig",
        qlow = longest_contig$range[1],
        qhigh = longest_contig$range[2]
      ),
      N50 = data.frame(
        name = "N50",
        qlow = N50$range[1],
        qhigh = N50$range[2]
      ),
      N95 = data.frame(
        name = "N95",
        qlow = N95$range[1],
        qhigh = N95$range[2]
      ),
      n_count = data.frame(
        name = "n_count",
        qlow = n_count$range[1],
        qhigh = n_count$range[2]
      ),
      kraken2_taxid_freq = data.frame(
        name = "kraken2_taxid_freq",
        qlow = kraken2_taxid_freq$range[1],
        qhigh = kraken2_taxid_freq$range[2]
      )
    )
    
    print(thresholds)
    
    # gc_content
    
    print("Keep rows where GC content is not too high or too low")
    
    if (length(which(qc$gc_content == FALSE)) > 0) {
      df_rm <- data.frame(
        assembly = qc$assembly[which(qc$gc_content == FALSE)],
        reason = "GC content is too high or too low"
      )
      df_failed <- dplyr::bind_rows(
        df_failed,
        df_rm
      )
      df$qc_pass[which(df$assembly %in% df_rm$assembly)] <- FALSE
    }
    
    print("QC passed overall:")
    print(sum(df$qc_pass))
    
    # contig_count
    
    print("Keep rows where contig count is not too high")
    
    if (length(which(qc$contig_count == FALSE)) > 0) {
      df_rm <- data.frame(
        assembly = qc$assembly[which(qc$contig_count == FALSE)],
        reason = "Contig count is too high"
      )
      df_failed <- dplyr::bind_rows(
        df_failed,
        df_rm
      )
      df$qc_pass[which(df$assembly %in% df_rm$assembly)] <- FALSE
    }
    
    print("QC passed overall:")
    print(sum(df$qc_pass))
    
    # longest_contig
    
    print("Keep rows where longest contig is not too short")
    
    if (length(which(qc$longest_contig == FALSE)) > 0) {
      df_rm <- data.frame(
        assembly = qc$assembly[which(qc$longest_contig == FALSE)],
        reason = "Longest contig is too short"
      )
      df_failed <- dplyr::bind_rows(
        df_failed,
        df_rm
      )
      df$qc_pass[which(df$assembly %in% df_rm$assembly)] <- FALSE
    }
    
    print("QC passed overall:")
    print(sum(df$qc_pass))
    
    # N50
    
    print("Keep rows where N50 is not too short")
    
    if (length(which(qc$N50 == FALSE)) > 0) {
      df_rm <- data.frame(
        assembly = qc$assembly[which(qc$N50 == FALSE)],
        reason = "N50 is too short"
      )
      df_failed <- dplyr::bind_rows(
        df_failed,
        df_rm
      )
      df$qc_pass[which(df$assembly %in% df_rm$assembly)] <- FALSE
    }
    
    print("QC passed overall:")
    print(sum(df$qc_pass))
    
    # N95
    
    print("Keep rows where N95 is not too short")
    
    if (length(which(qc$N95 == FALSE)) > 0) {
      df_rm <- data.frame(
        assembly = qc$assembly[which(qc$N95 == FALSE)],
        reason = "N95 is too short"
      )
      df_failed <- dplyr::bind_rows(
        df_failed,
        df_rm
      )
      df$qc_pass[which(df$assembly %in% df_rm$assembly)] <- FALSE
    }
    
    print("QC passed overall:")
    print(sum(df$qc_pass))
    
    # n_count
    
    print("Keep rows where N count is not too high")
    
    if (length(which(qc$n_count == FALSE)) > 0) {
      df_rm <- data.frame(
        assembly = qc$assembly[which(qc$n_count == FALSE)],
        reason = "N count is too high"
      )
      df_failed <- dplyr::bind_rows(
        df_failed,
        df_rm
      )
      df$qc_pass[which(df$assembly %in% df_rm$assembly)] <- FALSE
    }
    
    print("QC passed overall:")
    print(sum(df$qc_pass))
    
    # kraken2_taxid_freq
    
    print("Keep rows where kraken2_taxid_freq is not too low")
    
    if (length(which(qc$kraken2_taxid_freq == FALSE)) > 0) {
      df_rm <- data.frame(
        assembly = qc$assembly[which(qc$kraken2_taxid_freq == FALSE)],
        reason = "kraken2_taxid_freq is too low"
      )
      df_failed <- dplyr::bind_rows(
        df_failed,
        df_rm
      )
      df$qc_pass[which(df$assembly %in% df_rm$assembly)] <- FALSE
    }
    
    print("QC passed overall:")
    print(sum(df$qc_pass))
    
    ip1 <- sum(df$qc_pass)
    counter <- counter + 1
    
  }
  
  df_failed <- df_failed[order(df_failed$assembly), ]
  
  return(list(
    data = df,
    failed = df_failed,
    thresholds = thresholds
  ))
  
}
