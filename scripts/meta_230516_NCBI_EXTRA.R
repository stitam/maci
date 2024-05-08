# THIS IS AN UNTRACTED MODIFIED VERSION OF THE PIPELINE SCRIPT WITH SAME NAME.
# THIS SCRIPT CREATES THE METADATA TABLE FOR ASSEMBLIES THAT WERE DECLARED AS
# MISSING IN "220912_NCBI".

# This script requires two files downloaded from NCBI
# - assembly_result.xml: search for Acinetobacter baumannii within the NCBI
#   Assembly database:
#   https://www.ncbi.nlm.nih.gov/assembly/?term=acinetobacter%20baumannii
#   upper right corner -> send to -> file -> format = xml -> create file
# - biosample_result.txt: search for Acinetobacter baumannii within the NCBI
#   BioSample database:
#   https://www.ncbi.nlm.nih.gov/biosample/?term=acinetobacter+baumannii
#   upper right corner -> send to -> file -> format = full (text) -> create file
# - ncbi_strains_in_lab.xlsx
#
# This script requires the following extra files:
# - ncbi_strains_in_lab.tsv: a table which link internal ids and external ids
# - missing_assemblies_first_run.tsv: a table which lists genomes that should 
#   have assembly metadata and biosample metadata but the assembly file cannot
#   be found. These are caprured during input validation in main.nf.
# - missing_assemblies_second_run.tsv: a table which lists genomes that should 
#   have assembly metadata and biosample metadata but the assembly file cannot
#   be found. These are caprured during input validation in main.nf.
#
# This script generates a number of new files
# - geographic_locations.xlsx: see /data/geographic_locations.xlsx
# - geographic_locations_new.xlsx (optional)
# - geocodes.rda: see /data/geocodes.rda
# - non_human_related_isolates.tsv
# - metadata.rds
# - metadata.tsv

rm(list=ls())

##### SET JOB NAME #####

jobname <- "230516_NCBI_EXTRA"

########################

library(ggmap)
library(taxizedb)
library(testthat)

# TODO replace with library() once seqdb is public
library(devtools)
load_all("~/Methods/aci")
load_all("~/Methods/seqdb")

# parse assembly_result.xml into a tibble and export
if (!file.exists("assembly_meta.rds")) {
  assembly_meta <- ncbi_parse_assembly_xml("assembly_result.xml")
  
  # ncbi_parse_assembly_xml returns coverage as 'character' on purpose
  # because some values are '>10000'. However, validation requires this variable
  # to be a numeric.
  
  assembly_meta$coverage <- gsub(">", "", assembly_meta$coverage)
  assembly_meta$coverage <- as.numeric(assembly_meta$coverage)
  saveRDS(assembly_meta, file = "assembly_meta.rds")
} else {
  assembly_meta <- readRDS("assembly_meta.rds")
}

# parse biosample_result.txt into a tibble and export
if (!file.exists("biosample_meta.rds")) {
  biosample_meta <- ncbi_parse_biosample_txt(
    "biosample_result.txt", verbose = TRUE)
  saveRDS(biosample_meta, file = "biosample_meta.rds")
} else {
  biosample_meta <- readRDS("biosample_meta.rds")
}

# identify human related samples

# look at variable "attr_host"

include <- c(
  "Homo sapiens",
  "homosapiens",
  "hosptial surface", # intentional typo
  "Homosapiens",
  "Hemoculture",
  "Catheter",
  "Cerebrospinal fluid",
  "Humans",
  "sputum"
)

exclude <- c(
  "Murine",
  "Mus musculus",
  "Pediculus humanus humanus",
  "Lizard",
  "Environment",
  "Bos taurus",
  "Canis lupus familiaris",
  "Sus scrofa domesticus",
  "chicken",
  "environmental",
  "white stork",
  "white stork nestling",
  "Snake",
  "pig",
  "Apis mellifera",
  "dog",
  "duck",
  "goose",
  "natural / free-living",
  "cat",
  "Duck",
  "Sus scrofa",
  "Agama",
  "Canis lupus",
  "Cattle",
  "Dendrocygna viduata",
  "domestic grey parrot",
  "grey parrot",
  "Healthy coral",
  "horse",
  "Lettuce",
  "Mice",
  "Musca domestica",
  "Ovis aries",
  "Parthenium argentatum Gray (guayule shrubs)",
  "soil",
  "turkey",
  "Laboratory"
)

biosample_meta$human_related <- in_group(
  biosample_meta$attr_host,
  include, 
  exclude
)

# attr_geographic_location

include <- c(
  "China: Hangzhou, Sir Run Run Shaw Hospital ICU",
  "USA: Pittsburgh, Magee-Womens Hospital of UPMC",
  "USA: Strain isolated from hospital outbreak in Los Angeles County and provided by the Los Angeles County Public Health Laboratory."
)

exclude <- c(
  "Brazil: Lake Agua Preta, PA",
  "China: Huaishudian Metro Station, Chengdu",
  "Dominican Republic:Santo Domingo, Isabela River",
  "Russia: Tyumen region, river Iryum"
)

biosample_meta$human_related <- in_group(
  biosample_meta$attr_geographic_location,
  include,
  exclude,
  cache = biosample_meta$human_related
)

# "attr_isolate"

include <- c(
  "clinical"
)

exclude <- c(
  "Dog",
  "Cat",
  "Environmental isolate"
)

biosample_meta$human_related <- in_group(
  biosample_meta$attr_isolate,
  include,
  exclude,
  cache = biosample_meta$human_related
)

# "isolation_source"

include <- c(
  "wound",
  "Surveillance", # collected by United States Army, within the U.S. military healthcare system
  "human",
  "clinical isolates",
  "blood",
  "nosocomial",
  "sputum",
  "bedside rail in hospital intensive care unit",
  "respiratory",
  "urine",
  "Sterile Body Site",
  "alcohol foam dispenser in hospital intensive care unit",
  "bedside light switch in hospital intensive care unit",
  "nursing call button in hospital intensive care unit",
  "Tissue",
  "washroom sink in hospital intensive care unit",
  "hospital surface",
  "Warwound",
  "clinical isolate",
  "Fatal meningitis of a 4-month old infant",
  "Sterile Fluid",
  "Burned skin",
  "cerebrospinal fluid",
  "hospital bed",
  "Human urine",
  "Nosocomial",
  "Bronchial aspirate",
  "de-identified patient sample",
  "hospital drawer",
  "Hospital Effluent",
  "Medical waste water",
  "tracheal tube",
  "bronchial secretion",
  "catheter",
  "hospital",
  "Hospital Environment",
  "Hospital Floor",
  "Hospital ventilator",
  "hospital wastewater",
  "Human abscess",
  "bronchial aspirate",
  "bronchoalveolar lavage fluid",
  "Chair in hospital room",
  "Clinical sample, blood",
  "Clinical sample, endotracheal aspirate",
  "Clinical sample, sputum.",
  "clinical strain",
  "endotracheal aspiration",
  "eye after gunshot wound",
  "Fatal meningitis",
  "hospital chair",
  "Hospital effluent",
  "Hospital room surface",
  "hospital sewage",
  "hospital wastewaters",
  "human gut",
  "human skin",
  "Infected Wound",
  "isolated from the sputum of an intensive care unit patient at Meizhou People's Hospital",
  "Mucopurulent Sputum",
  "nasal swab (community carriage sample, isolated within 48 hours of hospital admission), Hospital for Tropical Diseases (Ho Chi Minh City)",
  "Phlegmone",
  "Punctate abdominal cavity",
  "Sewage sample from Ruijin Hospital",
  "Spinal fluid",
  "tracheal swab",
  "urinary",
  "Urinary culture",
  "urinary tract infection",
  "Urine sample sampled in I.R.C.C.S. Fondazione S. Maugeri neurological rehabilitation unit",
  "wound care scissors",
  "Wound swab",
  "bedside table", # collected by Clinical Hospital from Federal University, Brazil
  "Room 7", # collected by Concord Repatriation Hospital Sydney Australia
  "bloob", # isolate is human
  "Cavum" # collected by Nini Hospital, Lebanon
)

exclude <- c(
  "Laboratory evolved biofilm derived from ATCC17978", 
  "Draining-Matting Conveyor",
  "Laboratory evolved culture derived from ATCC17978",
  "SALMON (SA)",
  "wastewater",
  "ENVIRONMENT",
  "SHRIMP (SH)",
  "laboratory",
  "Lab",
  "environmental isolate",
  "ground",
  "Isabela River",
  "laboratory culture",
  "sewage",
  "soil",
  "pets",
  "school toilet",
  "Sewage",
  "water", 
  "food",
  "waste water",
  "wood",
  "raw pet food",
  "Animal-Calf-Bob Veal",
  "Environment",
  "field soil",
  "Manis javanica",
  "mouse",
  "Ozama River",
  "paper mill kaolin",
  "Raw sewage",
  "river water",
  "school toilet floor",
  "a feces sample of chicken origin",
  "activated sludge",
  "Activated sludge from the top of an aeration basin at the Des Moines Metropolitan Wastewater Reclamation Authority",
  "Air conditioner veterinary clinic",
  "Andrias davidianus",
  "Animal (cat)",
  "animal-cattle-steer",
  "cucumber rhizosphere",
  "Enviromental",
  "environment",
  "environmental",
  "Environmental sample pools",
  "inner root",
  "mangrove soil",
  "mouse gut",
  "Mustard waste water",
  "oil reservoir",
  "paper mill pulp",
  "paper pulp mill",
  "plant microbiome",
  "Plant microbiota",
  "river",
  "River Bhagirathi",
  "roots from fresh vegetable",
  "Seine River",
  "Shenzhou11 spacecraft",
  "Shivkund Water sample",
  "soil sample",
  "Tankmilk",
  "Tiangong2 space laboratory",
  "vegetable market soil",
  "Wastewater influent sample",
  "wastewater treatment plants",
  "Water (artesian well)",
  "water sample around Tzu Chi University, Hualien, Taiwan",
  "lab derived",
  "laboratory-derived",
  "metal/plastic", # metagenomic assembly obtained from the terrestrial metagenome reads
  "plastic", # terrestrial metagenome
  "kitchen counter", # metagenomic assembly
  "metal", # metagenomic assembly
  "In vitro experiment", # Polyamine shock RNAseq sample spermine treatment, Australia
  "Luria Broth Agar containing Benzalkonium Chloride", # BZK mutants, USA
  "Exponentially Growing Culture in LB broth", # WT, gigA, gigB mutants grown in lab, USA
  "Kanamycin Exposed Culutre in LB broth", # WT, gigA, gigB mutants grown in lab, USA
  "FDA-Stanford-Microbiome isolate"
)

biosample_meta$human_related <- in_group(
  biosample_meta$attr_isolation_source,
  include,
  exclude,
  cache = biosample_meta$human_related
)

# attr_isolate_souce

include <- "hemoculture"
exclude <- character()

biosample_meta$human_related <- in_group(
  biosample_meta$attr_isolate_souce,
  include,
  exclude,
  cache = biosample_meta$human_related
)

# attr_collected_by

include <- c(
  "CDC", # Pathogen: clinical or host-associated sample from Acinetobacter baumannii
  "Clinical Microbiology Department",
  "United States Army",
  "Nini Hospital"
)

exclude <- c(
  "FDA-Stanford Microbiome Effort",
  "Veterinary Research Institute"
)

biosample_meta$human_related <- in_group(
  biosample_meta$attr_collected_by,
  include,
  exclude,
  cache = biosample_meta$human_related
)

# broad.scale_environmental_context

include <- c(
  "human skin",
  "tissue",
  "Human",
  "human-associated",
  "human"
)

exclude <- c(
  "Surface sea water",
  "wastewater/sludge",
  "Lab culture",
  "synthetic metagenome"
)

biosample_meta$human_related <- in_group(
  biosample_meta$attr_broad.scale_environmental_context,
  include,
  exclude,
  cache = biosample_meta$human_related
)

# attr_insdc_center_name

include <- "University Hospital Muenster"
exclude <- character()

biosample_meta$human_related <- in_group(
  biosample_meta$attr_insdc_center_name,
  include,
  exclude,
  cache = biosample_meta$human_related
)

# title

include <- c(
  "Pathogen: clinical or host-associated sample from Acinetobacter baumannii",
  "ESBL-producing Enterobacteriaceae from neonates in Tanzania"
)

exclude <- c(
  "swab with mock bacterial cells, incubated at temperature",
  "Antibiotic-resistant bacterial isolates from water",
  "Acinetobacter isolates from raw meat in Singapore",
  "Metagenome or environmental sample from synthetic metagenome",
  "Pathogen: environmental/food/other sample from Acinetobacter baumannii"
)

biosample_meta$human_related <- in_group(
  biosample_meta$title,
  include,
  exclude,
  cache = biosample_meta$human_related
)

# identify known outbreaks

biosample_meta$outbreak <- grepl("[O|o]utbreak", biosample_meta$title) |
  grepl("[O|o]utbreak", biosample_meta$description)

# filter to relevant variables
meta <- biosample_meta[,which(names(biosample_meta) %in% c(
  "biosample",
  "attr_collection_date",
  "attr_geographic_location",
  "attr_latitude_and_longitude",
  "attr_strain",
  "attr_host",
  "human_related",
  "outbreak"
))]

meta <- dplyr::rename(
  meta, 
  collection_date = attr_collection_date,
  geographic_location = attr_geographic_location,
  latitude_and_longitude = attr_latitude_and_longitude,
  strain = attr_strain,
  host = attr_host
)

# merge assembly and biosample metadata
meta <- dplyr::left_join(assembly_meta, meta)

# convert NA strings to NA
na_strings <- c(
  "",
  "not available",
  "Not available"
)

meta$collection_date <- ifelse(
  meta$collection_date %in% na_strings,
  NA,
  meta$collection_date
)

meta$jobname <- jobname

meta$assembly_source <- "ncbi"

# Export a table of geographic locations, edit manually
# /data/geographic_locations.xlsx contains the table for the downloaded genomes.
t <- data.frame(geographic_location = sort(unique(meta$geographic_location)))
t$geographic_location <- trimws(t$geographic_location)
if(!file.exists("geographic_locations.xlsx")){
  writexl::write_xlsx(t, path = "geographic_locations.xlsx")
  stop("Exported 'geographic_locations.xlsx'. Resolve manually.")
} else {
  geoloc <- readxl::read_excel(
    "geographic_locations.xlsx", na = "NA")
  index <- which(
    t$geographic_location %in% geoloc$geographic_location == FALSE)
  if (length(index) > 0) {
    geoloc <- dplyr::full_join(geoloc, t, by = "geographic_location")
    geoloc <- geoloc[order(geoloc$geographic_location),]
    # If there are new geographic locations, export the new locations so they
    # can be added to the original table manually.
    writexl::write_xlsx(geoloc, path = "geographic_locations_new.xlsx")
    stop("Exported 'geographic_locations_new.xlsx'. Resolve manually.")
  }
}

# parse geographic locations to search terms for ggmap

geoloc$query <- NA
for(i in 1:nrow(geoloc)) {
  index <- which(!is.na(geoloc[i, 1:5]))[-1]
  if (length(index) > 0) {
    index_rev <- sort(index, decreasing = TRUE)
    geoloc$query[i] <- paste(geoloc[i, index_rev], collapse = ", ")
  }
}
geoloc$query <- ifelse(
  geoloc$query == "alaska, usa",
  "alaska",
  geoloc$query
)

searchterms <- unique(geoloc$query)
searchterms <- searchterms[which(!is.na(searchterms))]

# query GPS coordinates using the Google Maps API
# /data/geocodes.rds contains the gps coordinates for the downloaded genomes
gps <- query_geocode(searchterms, cache_file = "geocodes")

geoloc <- dplyr::left_join(
  geoloc,
  gps,
  by = "query"
)

geoloc <- geoloc[, which(names(geoloc) %in% c(
  "geographic_location",
  "country",
  "city",
  "lon",
  "lat"
))]

meta$geographic_location <- trimws(meta$geographic_location)
meta <- dplyr::left_join(meta, geoloc, by = "geographic_location")

# Query ISO2C country codes
meta$country_iso2c <- countrycode::countrycode(
  meta$country, "country.name", c("iso2c", "ecb"))
meta$country_iso2c[which(meta$country == "south_korea")] <- "KR"

# Query continents
meta$continent <- tolower(countrycode::countrycode(
  meta$country, origin = "country.name", destination = "continent"))
meta$continent <- ifelse(
  meta$country == "south_korea", "asia", meta$continent)
meta$continent <- ifelse(
  meta$country == "kosovo", "europe", meta$continent)

# Query regions
meta$region23 <- countrycode::countrycode(meta$country, "country.name", "region23")
meta$region23 <- tolower(gsub(" +", "_", meta$region23))
meta$region23[which(meta$country == "kosovo")] <- "southern_europe"
meta$region23[which(meta$country == "south_korea")] <- "eastern_asia"

genomes <- list.files("../assembled_genomes")
genomes <- gsub("\\.fna", "", genomes)

# 2022-10-26 CHANGE
# Assemblies without metadata are no longer included in analysis

# index <- which(genomes %in% meta$assembly == FALSE)

# missing <- data.frame(
#   assembly = genomes[index],
#   jobname = spath[length(spath)]
# )

# meta <- dplyr::bind_rows(meta, missing)

# query organism (species) name from taxids

taxids <- unique(meta$taxid)
taxids <- taxids[which(!is.na(taxids))]

# taxizedb may require that the ncbi taxonomy database is downloaded firs
# taxizedb::db_download_ncbi()

taxa <- taxizedb::classification(taxids)

species <- data.frame(
  taxid = names(taxa),
  organism = sapply(taxa, function(x) x$name[which(x$rank == "species")])
)

meta <- dplyr::left_join(meta, species, by = "taxid")

meta <- meta[which(meta$organism == "Acinetobacter baumannii"),]

meta <- dplyr::distinct(meta)

meta$relative_path <- paste0(
  "assembled_genomes/", meta$assembly, ".fna")

# GCA/GCF filtering
# TODO rename function in seqdb and update name here
meta <- meta[which(seqdb::flag_files(meta$assembly)$keep == TRUE),]

# RENAME NCBI GENOMES WHICH WE ALSO HAVE IN THE LAB
# TODO this should be performed after merging analysis results with metadata
in_lab <- read.csv("ncbi_strains_in_lab.tsv", sep = "\t")

# Check that all external ids are in the database
for (i in in_lab$external_id) {
  expect_true(i %in% meta$assembly)
}

# ASSEMBLY FILE SHOULD BE ALREADY RENAMED IN CENTRAL REPOSITORY
# THE CODE BELOW ONLY RENAMES meta$assembly and meta$relative_path

# TODO check
meta$in_lab <- FALSE
for (i in 1:nrow(in_lab)) {
  index <- which(meta$strain == in_lab$strain[i])
  meta$assembly[index] <- paste0(
    in_lab$internal_id[i], "_", meta$assembly[index])
  meta$relative_path[index] <- paste0(
    "assembled_genomes/", meta$assembly[index], ".fna")
  meta$in_lab[index] <- TRUE
}

# metadata may contain assemblies that cannot be found within the downloads
# input validation for the 'aci' pipeline will return such a list
# only keep entries where both metadata and assembly are available
missing1 <- read.csv("missing_assemblies_first_run.tsv", sep = "\t")
missing2 <- read.csv("missing_assemblies_second_run.tsv", sep = "\t")
missing <- rbind(missing1, missing2)

##### THIS IS THE UNTRACKED MODIFICATION #####


# INSTEAD OF REMOVING THESE ROWS, ONLY KEEP THESE
# meta <- meta[which(meta$assembly %in% missing$assembly == FALSE),]

# EDIT 2023-05-18
# Add GCF_021569215.1 this is an EURECA sample which was not downloaded for some reason
meta <- meta[which(meta$assembly %in% c(missing$assembly, "GCF_021569215.1") == TRUE),]

# rename entries that were queried as GCA but GCF was found

ids <- c(
  "GCA_025791815.1",
  "GCA_025791795.1",
  "GCA_025792335.1",
  "GCA_025792315.1",
  "GCA_025791875.1",
  "GCA_025818195.1",
  "GCA_025818205.1",
  "GCA_025818235.1",
  "GCA_025791855.1",
  "GCA_025826995.1"
)

for (i in ids) {
  index <- which(meta$assembly == i)
  meta$assembly[index] <- gsub("GCA", "GCF", meta$assembly[index])
  meta$relative_path[index] <- gsub("GCA", "GCF", meta$relative_path[index])
}

# remove entries that could not be downloaded from NCBI

ids <- c(
  "GCA_947045005.1", 
  "GCA_947045045.1", 
  "GCA_025758815.1", 
  "GCA_025758775.1", 
  "GCA_025758755.1", 
  "GCF_000248255.1"
)

index <- which(meta$assembly %in% ids)

if (length(index) > 0) {
  meta <- meta[-index, ]
}

###################################

# remove files which input validation found too small

# small1 <- read.csv("smallfiles_first_run.tsv", sep = "\t")
# meta <- meta[which(meta$assembly %in% small1$assembly == FALSE),]

#meta <- meta[-which(meta$assembly %in% c(
#  "GCA_900519095.1", "GCA_900496275.1", "GCA_900496215.1", "GCA_900496035.1"
#)),]

meta <- meta[, -which(names(meta) %in% c(
  "asm_name",
  "bioproject",
  "wgs_project",
  "status",
  "taxid",
  "host"
))]

meta_val <- validate_metadata(meta, keep_extra = TRUE)

biosample_non_human <- meta_val$biosample[which(meta_val$human_related == FALSE)]

biosample_non_human <- biosample_meta[which(biosample_meta$biosample %in% biosample_non_human),]

entry_count <- apply(biosample_non_human, 2, function(x) sum(!is.na(x)))

biosample_non_human <- biosample_non_human[,which(entry_count > 0)]

write.table(
  biosample_non_human,
  file = "non_human_related_isolates.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

saveRDS(meta_val, file = "metadata.rds")

write.table(meta_val,
            file = "metadata.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
