library(devtools)
library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to pipeline directory."
  ),
  make_option(
    c("-a", "--aci_path"),
    type = "character",
    help = "Path to aci prediction results."
  ),
  make_option(
   c("-o", "--pp_path"),
   type = "character",
   help = "Path to PopPUNK prediction results."
  ),
  make_option(
    c("-c", "--collapse_strategy"),
    type = "character",
    help = "Collapse strategy, e.g. geodate."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "~/Methods/aci",
    aci_path = "aci_filtered.rds",
    pp_path = "ppdb_clusters.csv",
    collapse_strategy = "geodate"
  )
}

# import custom functions
load_all(args$project_dir)

# import data set
aci <- readRDS(args$aci_path)

if (args$collapse_strategy == "none") {
  saveRDS(aci, file = "aci_collapse_none.rds")
}

set.seed(0)

if (args$collapse_strategy == "geodate") {
  
  ##### Option 1 - collapse by location and date #####
  
  # Strategy for collapsing outbreaks:
  
  # If isolates from the same city and the same year and month have the same ST-K
  # combined serotype, then they are considered to be samples from the same
  # outbreak and only one of them is kept per month, at random.
  
  # For isolates where the city is not known, if isolates from the same country 
  # and from the same year and same week of the year have the same ST-K combined
  # serotype, then they are considered to be samples from the same outbreak and
  # only one of them is kept per week, at random.
  
  # Isolates with unknown collection dates will not be eliminated from the data
  # set.

  # add new temporal variables
  aci$yearweek <- paste(
    aci$collection_year, lubridate::week(aci$collection_day), sep = "_")
  aci$yearweek <- ifelse(grepl("NA", aci$yearweek), NA, aci$yearweek)
  
  aci$yearmonth <- paste(
    aci$collection_year, lubridate::month(aci$collection_day), sep = "_")
  aci$yearmonth <- ifelse(grepl("NA", aci$yearmonth), NA, aci$yearmonth)
  
  # collapse samples where city is known
  cities <- unique(aci$city)[which(!is.na(unique(aci$city)))]
  keep_city <- vector()
  for (i in cities) {
    geo <- aci[which(aci$city == i),]
    serotypes <- unique(geo$serotype)
    for (j in serotypes) {
      geo_sero <- geo[which(geo$serotype == j),]
      for (k in unique(geo_sero$yearmonth)) {
        if (is.na(k)) {
          # if either year or month or both are not known, keep assembly
          keep_city <- c(
            keep_city, geo_sero$assembly[which(is.na(geo_sero$yearmonth))])
        } else {
          # otherwise keep one assembly per year per month
          geo_sero_yearmonth <- geo_sero[which(geo_sero$yearmonth == k),]
          keep_city <- c(keep_city, sample(geo_sero_yearmonth$assembly, 1))
        }
      }
    }
  }
  
  # collapse samples where city is not known but country is known
  countries <- unique(aci$country)[which(!is.na(unique(aci$country)))]
  keep_country <- vector()
  for (i in countries) {
    # only consider samples where the country is known but the city is not
    geo <- aci[which(aci$country == i & is.na(aci$city)),]
    serotypes <- unique(geo$serotype)
    for (j in serotypes) {
      geo_sero <- geo[which(geo$serotype == j),]
      for (k in unique(geo_sero$yearweek)) {
        if (is.na(k)) {
          # if either year or week or both are not known, keep assembly
          keep_country <- c(
            keep_country, geo_sero$assembly[which(is.na(geo_sero$yearweek))])
        } else {
          # otherwise keep one assembly per year per week
          geo_sero_yearweek <- geo_sero[which(geo_sero$yearweek == k),]
          keep_country <- c(keep_country, sample(geo_sero_yearweek$assembly, 1))
        }
      }
    }
  }
  
  index  <- which(
    aci$assembly %in% keep_city | 
      aci$assembly %in% keep_country |
      is.na(aci$country)
  )
  
  aci_collapsed <- aci[index, ]
  
  if (!interactive()) {
    saveRDS(aci_collapsed, file = "aci_collapse_geodate.rds")
  }
}

if (args$collapse_strategy == "poppunk") {
  
  ##### Option 2 - collapse by PopPUNK clusters #####
  
  # Strategy fo collapsing outbreaks:
  
  # Isolates are assigned to a PopPUNK cluster. The core accessory distance
  # thresholds are adjusted such that only very closely related isolates will be
  # assigned to the same PopPUNK cluster. Isolates that belong to the same popPUNK
  # cluster are considered to be from the same outbreak and only one of them is
  # kept per year.
  
  # TODO: Wondering if this is too strict. When isolates from the same year are
  # assigned to the same PopPUNK cluster, only one of them will be kept regardless
  # of their geographical location. If a cluster is present in multiple cities or
  # countries, this variability will be lost during downsampling.
  
  # add poppunk clusters
  pp <- read.csv(args$pp_path)
  pp <- dplyr::rename(pp, assembly = Taxon, pp = Cluster)
  pp$assembly <- gsub("_", ".", pp$assembly)
  pp$assembly <- gsub("GCA\\.", "GCA_", pp$assembly)
  pp$assembly <- gsub("GCF\\.", "GCF_", pp$assembly)
  aci <- dplyr::left_join(aci, pp, by = "assembly")
  
  # filter to genomes that have predicted poppunk clusters
  
  # TODO: Is this appropriate? The only way for an assembly not to have a PopPUNK
  # cluster assignment is that it was not included in the PopPUNK analysis. This
  # can happen if the assembly was not available (new sample) or if the assembly
  # was during PopPUNK quality control. Maybe samples that were dropped during
  # PopPUNK quality control should be dropped at the very beginning of the
  # analysis?
  
  aci <- aci[which(!is.na(aci$pp)), ]
  
  # keep one sample per poppunk cluster per country per year
  # ISSUE: poppunk clusters 1 and 2 contain multiple clusters, these will be significantly downsampled
  aci_ds <- data.frame()
  for (i in unique(aci$pp)) {
    aci_i <- aci[which(aci$pp == i),]
    for (j in unique(aci_i$country)) {
      aci_j <- aci_i[which(aci_i$country == j),]
      for (k in unique(aci_j$collection_year)) {
        aci_k <- aci_j[which(aci_j$collection_year == k),]
        aci_ds <- dplyr::bind_rows(aci_ds, aci_k[sample(1:nrow(aci_k), 1),])
      }
    }
  }
  
  if (!interactive()) {
    saveRDS(aci_ds, "aci_collapse_poppunk.rds")
  }
  
}
