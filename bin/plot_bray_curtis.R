library(optparse)
library(tidyverse)
rm(list = ls())

args_list <- list(
  make_option(
    "--project_dir",
    type = "character",
    help = "Path to project directory."
  ),
  make_option(
    c("-f", "--file"),
    type = "character",
    help = "Path to aci prediction results."
  ),
  make_option(
    c("-c", "--country_file"),
    type = "character",
    help = "Path to a tbl prepared for serotop by countries, e.g. top_serotypes_country_ds_geodate.tsv."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "aci",
    file = "aci_crab_ds_geodate.tsv",
    country_file = "top_serotypes_country_ds_geodate.tsv"
  )
}

library(devtools)
load_all(args$project_dir)

data <- read_delim(args$country_file)

# prepare the country_cont object from scratch
# country_cont <- read_csv("country_cont.csv")
aci <- read_df(args$file) %>% 
  dplyr::filter(filtered & crab & downsampled & downsampled_by_pop)

country_cont <- data.frame(
  continent = aci$continent,
  Var1 = aci$country,
  Var2 = aci$country,
  region23 = aci$region23
) %>% distinct()

country_serotypes_geodate <- read_delim(args$country_file)
# country_serotypes_geodate$country <- gsub("_"," ", country_serotypes_geodate$country)
# country_serotypes_geodate$country <- tools::toTitleCase(country_serotypes_geodate$country)

##### MANUAL SELECTION OF COUNTRIES #####
# rarefaction sample size > 40 in these European countries)
# TODO: this should either be done automatically or the list should be a
# pipeline argument
eu <- c("france", "serbia", "germany", "romania", "hungary", "italy", "greece")
# check if all countries are present in the data
if (any(!eu %in% country_serotypes_geodate$country)) {
  index <- which(!eu %in% country_serotypes_geodate$country)
  missing_collapsed <- paste(eu[index], collapse = ", ")
  msg <- paste0(
    "Not all countries in the EU list are present in the data: ",
    missing_collapsed
  )
  stop(msg)
}
#########################################

# between european countries
bray <- country_serotypes_geodate  %>%  
  filter(country_serotypes_geodate$country %in% eu) %>% 
  group_by(country, serotype) %>% 
  summarise(count=sum(count))

bray_wide <- pivot_wider(bray, names_from=serotype, values_from=count)
diversity_matrix <- as.matrix(bray_wide[,-1])
rownames(diversity_matrix) <- bray_wide$country
diversity_matrix[is.na(diversity_matrix)] <- 0
beta_diversity <- vegan::vegdist(diversity_matrix, "bray")

within_europe <- reshape2::melt(as.matrix(beta_diversity))
within_europe <- within_europe %>% filter(Var1!=Var2) %>% distinct(value, .keep_all = T) %>% subset(select=value) %>% mutate(type="within_europe")

# with non-European countries
bray2 <- country_serotypes_geodate %>% 
  group_by(country, serotype) %>% 
  summarise(count=sum(count))

bray2_wide <- pivot_wider(bray2, names_from=serotype, values_from=count)
diversity_matrix2 <- as.matrix(bray2_wide[,-1])
rownames(diversity_matrix2) <- bray2_wide$country
diversity_matrix2[is.na(diversity_matrix2)] <- 0
beta_diversity2 <- vegan::vegdist(diversity_matrix2, "bray")

btw_europe_and_other <- reshape2::melt(as.matrix(beta_diversity2))
btw_europe_and_other <- country_cont %>% subset(select=c(Var1, continent)) %>% left_join(btw_europe_and_other, country_cont, by="Var1")
colnames(btw_europe_and_other)[2] <- "continent1"
btw_europe_and_other <- country_cont %>% subset(select=c(Var2, continent)) %>% left_join(btw_europe_and_other, country_cont, by="Var2")
colnames(btw_europe_and_other)[2] <- "continent2"

btw_europe_and_other <- btw_europe_and_other %>% filter(continent2=="europe" & continent1!="europe")

##### MANUAL SELECTION OF COUNTRIES #####
# TODO: this should either be done automatically or the list should be a
# pipeline argument
non_eu <- c(
  "south_africa", "brazil", "thailand", "china", "usa", "saudi_arabia", "israel")
# check if all countries are present in the data
if (any(!non_eu %in% btw_europe_and_other$Var1)) {
  index <- which(!non_eu %in% country_serotypes_geodate$country)
  missing_collapsed <- paste(non_eu[index], collapse = ", ")
  msg <- paste0(
    "Not all countries in the non-EU list are present in the data: ",
    missing_collapsed
  )
  stop(msg)
}
#########################################

btw_europe_and_other <- btw_europe_and_other %>% filter(btw_europe_and_other$Var1 %in% non_eu & btw_europe_and_other$Var2 %in% eu) %>% 
  subset(select=value) %>% mutate(type="btw_europe_and_other")

# between regions
all <- c(eu, non_eu)

bray3 <- country_serotypes_geodate %>% 
  filter(country_serotypes_geodate$country %in% all) %>% 
  group_by(country, serotype) %>% 
  summarise(count=sum(count))

bray3_wide <- pivot_wider(bray3, names_from=serotype, values_from=count)
diversity_matrix3 <- as.matrix(bray3_wide[,-1])
rownames(diversity_matrix3) <- bray3_wide$country
diversity_matrix3[is.na(diversity_matrix3)] <- 0
beta_diversity3 <- vegan::vegdist(diversity_matrix3, "bray")

btw_region <- reshape2::melt(as.matrix(beta_diversity3))
btw_region <- country_cont %>% subset(select=c(Var1, region23)) %>% left_join(btw_region, country_cont, by="Var1")
colnames(btw_region)[2] <- "region1"
btw_region <- country_cont %>% subset(select=c(Var2, region23)) %>% left_join(btw_region, country_cont, by="Var2")
colnames(btw_region)[2] <- "region2"

btw_region <- btw_region %>% filter(region2!=region1) %>% distinct(value) %>% mutate(type="btw_region")

data <- rbind(within_europe, btw_europe_and_other, btw_region)
data <- transform(data, value=as.numeric(value))

data$type <- factor(data$type, levels=c("within_europe", "btw_europe_and_other", "btw_region"))

# set seed for ggplot2::position_jitter()
set.seed(0)

g <- ggplot(data, aes(type, value)) +
  geom_boxplot(outlier.color="transparent") +
  geom_point(position = position_jitter(0.1, seed = 0)) +
  ggpubr::geom_signif(
    comparisons=list(
      c("within_europe", "btw_europe_and_other"),
      c("within_europe", "btw_region"),
      c("btw_europe_and_other", "btw_region")
    ),
    map_signif_level=TRUE,
    step_increase = 0.1,
    test="wilcox.test"
  ) +
  xlab("") +
  ylab("Bray-Curtis dissimilarity index")

##### NOT USED #####
# TODO do we need this?
a <- wilcox.test(within_europe$value, btw_europe_and_other$value, paired = F)
b <- wilcox.test(within_europe$value, btw_region$value, paired = F)
c <- wilcox.test(btw_europe_and_other$value, btw_region$value, paired = F)
labels=c(
  signif(a$p.value, 4), 
  signif(b$p.value, 4),
  signif(c$p.value, 4)
)
###################

ggsave(
  filename = "bray_curtis.pdf",
  plot = g,
  width = 8,
  height = 6
)
ggsave(
  filename = "bray_curtis.png",
  plot = g,
  width = 8,
  height = 6
)
saveRDS(g, "bray_curtis.rds")