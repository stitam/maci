#install.packages("stringi")
library(tidyverse)
library(stringi)

rm(list=ls())

#world regions

tbl_r<- read_tsv("top_serotypes_region23_ds_geodate.tsv")

#tbl_r %>% distinct(region23)

tbl_r$region23 <- stri_replace_all_regex(tbl_r$region23, c('eastern_europe','southern_europe','southern_asia', 'northern_america','south-eastern_asia','eastern_asia', 'western_asia','western_europe','south_america'), 
                                        c('Eastern Europe','Southern Europe','Southern Asia', 'Northern America','South-eastern Asia','Eastern Asia', 'Western Asia','Western Europe','South America'), vectorize=F)


tbl_r90 <- tbl_r %>%
  group_by(region23) %>% 
  filter(ratio_cumsum<0.90) %>% summarise(n=n()+1) 

g <- tbl_r90 %>%
  ggplot(mapping = aes(x=region23, y=n))+
  geom_bar(stat = "identity", color = "darkgray", fill = "darkgray")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle=90, color="black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  scale_x_discrete(limits=c("Northern America","South America","Western Europe","Southern Europe","Eastern Europe","Western Asia", "Southern Asia", "South-eastern Asia", "Eastern Asia"))+
  xlab("") +
  ylab("Nr. of top serotypes")

ggsave(
  filename = "Nr_top90_serotypes_regions.pdf",
  plot = g,
  width = 8,
  height = 6
)
ggsave(
  filename = "Nr_top90_serotypes_regions.png",
  plot = g,
  width = 8,
  height = 6
)
saveRDS(g, "Nr_top90_serotypes_regions.rds")


# countries

rm(list=ls())


tbl_c<- read_tsv("top_serotypes_country_ds_geodate.tsv")

#tbl_c %>% distinct(country)

tbl_c$country <- stri_replace_all_regex(tbl_c$country, c('usa','brazil','germany','france','greece','italy','hungary','romania','serbia','israel','saudi_arabia', 'thailand','china','south_africa'), 
                                        c('USA','Brazil','Germany','France','Greece','Italy','Hungary','Romania','Serbia', 'Israel','Saudi Arabia', 'Thailand','China','South Africa'), vectorize=F)


tbl_c90 <- tbl_c %>%
  group_by(country) %>% 
  filter(ratio_cumsum<0.90) %>% summarise(n=n()+1) 

g <- tbl_c90 %>%
  ggplot(mapping = aes(x=country, y=n))+
  geom_bar(stat = "identity", color = "darkgray", fill = "darkgray")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle=90, color="black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  scale_x_discrete(limits=c("USA","Brazil","Germany","France","Greece","Italy","Hungary","Romania","Serbia", "Israel","Saudi Arabia", "Thailand","China","South Africa"))+
  xlab("") +
  ylab("Nr. of top serotypes")

ggsave(
  filename = "Nr_top90_serotypes_countries.pdf",
  plot = g,
  width = 8,
  height = 6
)
ggsave(
  filename = "Nr_top90_serotypes_countries.png",
  plot = g,
  width = 8,
  height = 6
)
saveRDS(g, "Nr_top90_serotypes_countries.rds")