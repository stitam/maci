# This script generates a world map with the number of CRAB and non-CRAB
# isolates per region. The script only includes regions with at least 10
# isolates.

library(dplyr)
library(ggplot2)
library(optparse)
library(tidyr)

rm(list = ls())

args_list <- list(
  make_option(
    c("-f", "--file"),
    type = "character",
    help = "Path to aci prediction results."
  ),
  make_option(
    c("-r", "--regions"),
    type = "character",
    help = "Path to a table with region names, coordinates, colors."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    file = "aci_filtered.rds",
    regions = "geographic_locations_in_study.tsv"
  )
}

# read prediction results
aci_all <- readRDS(args$file)
# get list of countries for consistency check
countries <- sort(unique(aci_all$country))
# get list of regions for consistency check
regions <- sort(unique(aci_all$region23))

# read regions database
geoloc <- read.csv(args$regions, sep = "\t", na = "")
geoloc_countries <- geoloc %>% filter(variable == "country")
geoloc_regions <- geoloc %>% 
  filter(variable == "region23") %>%
  dplyr::rename(region23 = term, region23_pretty = term_pretty)

# import world map
countries_map <- map_data("world")

# check all countries in database can be found on world map
if (any(geoloc_countries$term_maps %in% countries_map$region == FALSE)) {
  index_missing <- which(!geoloc_countries$term_maps %in% countries_map$region)
  missing_collapsed <- paste(geoloc_countries$term_maps[index_missing], collapse = ", ")
  msg <- paste0(
    "Some country names are missing from map database: ",
    missing_collapsed,
    ". Probably a typo in the geographic locations database."
  )
  stop(msg)
}

# Filter to regions with at least 10 isolates
regions <- table(aci_all$region23)
regions <- regions[which(regions >= 10)]
regions <- names(regions)

# Filter to isolates from selected regions
aci_all <- aci_all[which(aci_all$region23 %in% regions),]
countries <- sort(unique(aci_all$country))
geoloc_countries <- geoloc_countries %>% filter(term %in% countries)
geoloc_regions  <- geoloc_regions %>% filter(region23 %in% regions)

color_country <- subset(aci_all, select=c(region23, country)) %>% 
  dplyr::distinct() %>%
  left_join(
    .,
    subset(geoloc_regions, select = c(region23, colour)),
    by = "region23"
  ) %>% 
  subset(., select = c(country, colour)) %>%
  dplyr::rename(region = country)

color_country$region <- sapply(color_country$region, function(x) {
  index <- which(geoloc_countries$term == x)
  return(geoloc_countries$term_pretty[index])
})

region_lab_data <- countries_map %>%
  group_by(region) %>%
  summarise(long = mean(long), lat = mean(lat))

color_country <- left_join(region_lab_data, color_country, by="region")
color_country <- replace(color_country, is.na(color_country), "#e1e3e1")
countries_map2 <- left_join(countries_map, color_country, by="region")


# data for the charts
# aci_all2 <- aci_all[which(aci_all$region23 %in% aci$region23),]
chart <- aci_all %>% 
  group_by(region23, crab) %>% 
  summarise(count=length(unique(assembly))) %>%
  filter(region23 %in% regions)

chart <- left_join(chart, geoloc_regions, by= "region23")

chart$crab[chart$crab=="FALSE"] <- "Non-CRAB"
chart$crab[chart$crab=="TRUE"] <- "CRAB"

chart2 <- chart %>% subset(select=c(region23_pretty, crab, count, lon, lat)) %>% 
  pivot_wider(names_from=crab, values_from=count)

chart2$`Non-CRAB` <- ifelse(is.na(chart2$`Non-CRAB`), 0, chart2$`Non-CRAB`)
chart2$CRAB <- ifelse(is.na(chart2$CRAB), 0, chart2$CRAB)

chart2 <- mutate(chart2, total=CRAB+`Non-CRAB`)

# only keep regions with at least 10 isolates

chart2 <- chart2[which(chart2$total >= 10),]

# chart size adjustment
chart2$radius <- log2(chart2$total)-4

world_map <- map_data("world")

world_map_fig <- ggplot(world_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill=countries_map2$colour, colour = "gray40", linewidth=0.1)+
  theme_void()+
  scatterpie::geom_scatterpie(
    aes(x=lon, y=lat, group=region23_pretty, r=radius, fill=type),
    alpha=0.8,
    data = chart2, 
    cols=c("CRAB", "Non-CRAB"),
    linewidth = 0.1
  )+
  #scatterpie::geom_scatterpie_legend(chart2$radius, x=-130, y=-40, labeller = function(x) 2^(x+4))+
  scale_fill_manual(name="", values=c("seagreen4", "gray50"))+
  guides(fill = "none")+
  theme(
    legend.text = element_text(size=12),
    legend.position = c(0.193,0.45),
    plot.margin = unit(c(0,0,0,0), "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )


region_color <- geoloc_regions %>%
  subset(select = c(region23_pretty, colour))

lg <- grid::legendGrob(labels=region_color$region23_pretty, 
           pch=rep(15,length(region_color$region23_pretty)),
           gp=grid::gpar(col=region_color$colour, fill="gray"),
           byrow=F, hgap=0.5, vgap = 0.2, ncol = 8)

#g <- ggplotify::as.ggplot(gridExtra::grid.arrange(world_map_fig, lg,  nrow=2, heights=c(10, 1)))

g <- world_map_fig

g <- g + theme(
  legend.text = element_text(size = 5)
)

ggsave(
    filename = "Fig1A_world_map.pdf",
    plot = g,
    units = "mm",
    width = 60,
    height = 40
)

ggsave(
  filename = "Fig1A_world_map.png",
  plot = g,
  units = "mm",
  width = 60,
  height = 40
)

saveRDS(g, file = "Fig1A_world_map.rds")