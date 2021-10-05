rm(list=ls())

url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/",
              "Acinetobacter_baumannii/assembly_summary.txt")

parse_assembly_summary <- function(url) {
  t <- readLines(url)
  body <- lapply(t[3:length(t)], function(x) strsplit(x, split = "\t")[[1]])
  df <- as.data.frame(do.call(rbind, body)) 
  names(df) <- strsplit(gsub("# ", "", t[2]), split = "\t")[[1]]
  return(df)
}

ncbi <- parse_assembly_summary(url)
ncbi <- dplyr::rename(ncbi, strain = infraspecific_name)
ncbi$strain <- gsub("strain=", "", ncbi$strain)

local <- read.csv("../input/acinetobacter_serotypes_treesort.tsv",
                       header = TRUE, sep = "\t")


meta <- ncbi[which(ncbi$assembly_accession %in% local$Assembly),]

write.table(meta,
            file = "./output/assembly_summary.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

download_biosample <- function(biosample_id) {
  cmd <- paste0("efetch -db biosample -format full -id ",
                biosample_id," > ", "./biosamples/", biosample_id, ".txt")
  system(cmd)
}

#uncomment the line below to download biosample data again, will take a while
#sapply(meta$biosample, download_biosample)

parse_biosample <- function(file) {
  t <- readLines(file)
  
  if (sum(grepl("Error", t)) > 0) {
    df <- data.frame()
    return(df)
  }
  
  match <- t[grep("^ */.+=", t)]
  
  if(length(match) == 0) {
    df <- data.frame()
    return(df)
  }
  
  match <- gsub("/", "", match)
  match <- gsub('"', "", match)
  match <- trimws(match)
  
  names <- sapply(match, function(x) strsplit(x, split = "=")[[1]][1])
  names <- gsub(" ", "_", unname(names))
  values <- lapply(match, function(x) strsplit(x, split = "=")[[1]][2])
  values <- unname(values)
  attributes <- as.data.frame(values)
  names(attributes) = names
  
  df <- data.frame(
    biosample = strsplit(basename(file), split = ".txt")[[1]][1],
    attributes
  )
  return(df)
}

dirpath <- "./biosamples"
filepaths <- paste0(dirpath, "/", dir(dirpath))
newmeta_list <- lapply(filepaths, parse_biosample)
newmeta <- dplyr::bind_rows(newmeta_list)

write.table(newmeta,
            file = "./output/biosample_summary.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

newmeta <- dplyr::rename(newmeta, strain2 = strain)

meta2 <- plyr::join(meta, newmeta, by = "biosample", type = "left")

meta2$strain3 <- sapply(meta2$organism_name, function(x){
  strsplit(x, split = "Acinetobacter baumannii ")[[1]][2]
})

unique_strains <- function(x) {
  x <- x[x != ""]
  x <- x[!is.na(x)]
  x <- unique(sort(x))
  if (length(x) > 1) x <- paste(x, collapse = "|")
  return(x)
}

meta2$strains <- apply(meta2[,c("strain","strain2","strain3", "isolate")],
                       1,
                       unique_strains)

meta3 <- meta2[,c(
  "assembly_accession",
  "asm_name",
  "bioproject",
  "biosample",
  "strains",
  "geographic_location",
  "latitude_and_longitude",
  "collection_date",
  "host",
  "host_disease",
  "isolation_source"
)]

write.table(meta3,
            file = "./output/acinetobacter_metadata.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

meta3 <- read.csv("./output/acinetobacter_metadata.tsv",
                    sep = "\t", na.strings = c(
                      "",
                      "unknow","unknown","Unknown",
                      "missing", "Missing",
                      "not collected", "Not collected",
                      "not applicable", "Not applicable",
                      "not available", "Not available",
                      "not determined"))

meta3$collection_year <- stringi::stri_sub(meta3$collection_date,
                                           from = 1,
                                           to = 4)

meta3 <- meta3[,c(
  "assembly_accession",
  "asm_name",
  "bioproject",
  "biosample",
  "strains",
  "geographic_location",
  "latitude_and_longitude",
  "collection_date",
  "collection_year",
  "host",
  "host_disease",
  "isolation_source"
)]

meta <- meta3

save(meta, file = "./output/acinetobacter_metadata.rda")

write.table(meta3,
            file = "./output/acinetobacter_metadata.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
