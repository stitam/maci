library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-f", "--file"),
    type = "character",
    help = "Path to aci prediction results."
  ),
  make_option(
    c("-s", "--sensitivity"),
    type = "character",
    help = "Path to phage sensitivity lab results."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    file = "data/typing_summary_tables/aci_study.rds",
    sensitivity = "labdata/aci_project/Ab_all_strains_phages_spotassay_PFU.tsv"
  )
}

library(tidyverse)
library(scales)
library(ggpubr)
library(data.table)
library(vegan)
library(stringi)

###################################
#plotting based on ordered countries
###################################

aci <- readRDS(args$file)
aci$Strain <- unname(sapply(aci$assembly, function(x) {
  if (grepl("^GC[A|F]", x)) {
    out <- x
  } else {
    out <- strsplit(x, "_")[[1]][1]
  }
  return(out)
}))

aci_all <- read.csv(args$sensitivity, sep = "\t")
aci_all <- aci_all [, -1]

# merge without adding space

aci_all$Strain <- paste(aci_all$Name, aci_all$Identifier.of.the.strain, sep="")

# TODO: DANGER ZONE! DANGER ZONE! NAMES ARE SELECTED BY INDEX
phage_names <- names(aci_all)[7:ncol(aci_all)]

Aci_comb = left_join(x=aci_all,y=aci,by="Strain") 

# check that the number of rows has not been inflated by left_join()
testthat::expect_true(nrow(Aci_comb) == nrow(aci_all))

#Aci_comb %>% distinct(country)

# új oszlopok létrehozása 2 összevonásából

Aci_comb$strain_city <- paste(Aci_comb$Strain, "-", Aci_comb$city)

Aci_comb$ST_KL <- paste(Aci_comb$MLST, "-", Aci_comb$KL)

# Aci_comb %>% distinct(ST_KL)

# filter out those for which we don't know the ST_KL

Aci_comb <- Aci_comb %>% filter (ST_KL != "NA - NA")

# make the table shorter

ACI <- Aci_comb %>% select(c(
  "Strain",
  "Identifier.of.the.strain_long",
  "MLST",
  "KL",
  phage_names,
  "country",
  "city",
  "strain_city",
  "ST_KL"
))

###### filter for the 5 countries from the region we got samples and the most prevalent 10 serotypes plus the 16th

ctry3 <- ACI %>% filter(country == "hungary" | country == "romania" | country == "serbia" | country == "bosnia_and_herzegovina" | country == "montenegro") %>%
  filter (ST_KL == "ST2 - KL3"| ST_KL == "ST636 - KL40" | ST_KL == "ST492 - KL104" | ST_KL == "ST2 - KL2" |
            ST_KL == "ST1 - KL1" | ST_KL == "ST2 - KL9" | ST_KL == "ST2 - KL12" | ST_KL == "ST1 - KL17" |
            ST_KL == "ST2 - KL7" | ST_KL == "ST2 - KL32" | ST_KL == "ST2 - KL77")

# ctry3_F <- ctry3 %>% filter (!is.na (Highwayman) & !is.na (Silvergun) & !is.na (Fanak) & !is.na (PhT2-v2) & !is.na (Porter) & !is.na (Dino) & !is.na (Fishpie) & !is.na (Tama) & !is.na (Margaret) & !is.na (ABW132) &  
#!is.na (ABW311) & !is.na (Navy4-v2) & !is.na (Rocket) & !is.na (Konradin-v2) & !is.na (KissB))

# kiszűrni minden sort, ahol bármelyik fág NA, illetve fertőz, de nincs PFU (YES, Yes, yes) 
#összevonvan felírva függvénnyel

ctry3_F <- ctry3 %>% filter (if_all(c(Highwayman,Silvergun,Fanak,PhT2.v2,Porter,Dino,Fishpie,Tama,Margaret,ABW132,ABW311,Navy4.v2,Rocket,Konradin.v2,KissB), function(x) !is.na(x))) %>%
  filter (if_all(c(Highwayman,Silvergun,Fanak,PhT2.v2,Porter,Dino,Fishpie,Tama,Margaret,ABW132,ABW311,Navy4.v2,Rocket,Konradin.v2,KissB), function(x) !`%in%`(x, c("YES", "Yes", "yes"))))


ctry3_F$country <- stri_replace_all_regex(ctry3_F$country, c('hungary','romania','serbia', 'bosnia_and_herzegovina','montenegro'), 
                                          c('Hungary','Romania','Serbia', 'Bosnia and Herzegovina','Montenegro'), vectorize=F)

# filter out 6 strains from ST636 KL40 and 8 from ST2 KL12 to represent the strains more similar to their abundance (tbl_ST_KL)

ctry3_F <- ctry3_F %>% filter(Strain != "Aci337" & Strain != "Aci342" & Strain != "Aci127"  & Strain != "Aci200"  & Strain != "Aci361"  & Strain != "Aci362" ) %>% 
  filter(Strain != "Aci141" & Strain != "Aci143" & Strain != "Aci136"  & Strain != "Aci305"  & Strain != "Aci307"  & Strain != "Aci309" & Strain != "Aci260" & Strain != "Aci135" )

#tbl_ST_KL <- ctry3_F %>% 
  #group_by(ST_KL) %>%  
  #summarise(Nr = n()) %>%  
  #arrange(desc(Nr))

# TODO: DANGER ZONE! COLUMNS ARE SELECTED BY INDEX
ACI_long <- pivot_longer(data = ctry3_F, cols = 5:19, values_to = "value")


ACI_long$value[ACI_long$value=="NO"] <- "0"

ACI_long <- transform(ACI_long, value=as.numeric(value))

#ACI_long %>% distinct(country)
#ACI_long %>% distinct(ST_KL)


order_Strain <- ACI_long %>%
  mutate(country=factor(country, levels=c("Hungary", "Romania", "Serbia", "Bosnia and Herzegovina", "Montenegro" ))) %>% 
  arrange(country) %>%
  pull(Strain) %>% unique()

# prepare for the composit figure

library (patchwork)
library("viridis")
library("viridisLite")

no_of_colors <- 22
palette <- viridis_pal(option = "H")(no_of_colors)
palette
pie(rep(1, no_of_colors), col = palette, main = "viridis palette")

# I define 1, which was Lysis from without as 0

tidy_value <- function(x) {
  x <- as.numeric(x)
  x_floor <- floor(x)
  x_ceiling <- ceiling(x)
  if (x == 1) {
    return("0") 
  }
  if (x!= 0 & x != 1 & x <= 3) {
    stop("Detected low PFU that is currently not handled.")
  }
  if (x > 3 & x <= 6) {
    return("3-6")
  }
  if (x > 6) {
    return(paste0(x_ceiling - 1, "-", x_ceiling))
  }
  return(x)
}

ACI_long$value2 <- sapply(ACI_long$value, tidy_value)

# megadom a sorrendet, ami alapján a színeket beállítja
#ACI_long$value2 <- factor(ACI_long$value2, levels=c("0", "5-6","6-7","7-8", "8-9", "9-10", "10-11", "Lysis from without"))

fun_color_range <- colorRampPalette(c("cadetblue1", "steelblue")) # színskála létrehozása
my_colors <- fun_color_range(6) # 6-at használ a színskálából
my_colors2 <- c("white",my_colors)
names(my_colors2)<-c("0", "3-6","6-7","7-8", "8-9", "9-10", "10-11")

ACI_long$name <- factor(
  ACI_long$name,
  levels=c(
    "Highwayman", "Silvergun", "Fanak", "PhT2.v2", "Porter", "Dino", "Fishpie",
    "Tama", "Margaret" ,"ABW132", "ABW311", "Navy4.v2", "Rocket", "Konradin.v2",
    "KissB"
    ),
  labels = c(
    "Highwayman", "Silvergun", "Fanak", "PhT2-v2", "Porter", "Dino", "Fishpie",
    "Tama", "Margaret" ,"ABW132", "ABW311", "Navy4-v2", "Rocket", "Konradin-v2",
    "KissB"
    )
  )

ST_KL_names <- c(
  `ST2 - KL3` = "ST2\nKL3",
  `ST636 - KL40` = "ST636\nKL40",
  `ST492 - KL104` = "ST492\nKL104",
  `ST2 - KL2` = "ST2\nKL2",
  `ST1 - KL1` = "ST1\nKL1",
  `ST2 - KL9` = "ST2\nKL9",
  `ST2 - KL12` = "ST2\nKL12",
  `ST1 - KL17` = "ST1\nKL17",
  `ST2 - KL7` = "ST2\nKL7",
  `ST2 - KL32` = "ST2\nKL32",
  `ST2 - KL77` = "ST2\nKL77"
)

plot1 <- 
  ACI_long %>%
  mutate(name=factor(name, ordered = TRUE,levels=c("Highwayman", "Silvergun", "Fanak", "PhT2-v2", "Porter", "Dino", "Fishpie", "Tama", "Margaret" ,"ABW132",  
                                                   "ABW311", "Navy4-v2", "Rocket", "Konradin-v2", "KissB"))) %>% 
  mutate(value2=factor(value2, ordered = TRUE,levels = names(my_colors2)) )%>%
  mutate(Strain=factor(Strain, ordered = TRUE,levels = order_Strain )) %>% 
  mutate(ST_KL=factor(ST_KL, ordered = TRUE,levels=c("ST2 - KL3", "ST636 - KL40", "ST492 - KL104", "ST2 - KL2", "ST1 - KL1", "ST2 - KL9" , "ST2 - KL12" , "ST1 - KL17" ,
                                                     "ST2 - KL7" , "ST2 - KL32" , "ST2 - KL77"))) %>%
  ggplot(mapping = aes(x=Strain,y=name,  fill=value2))+
  geom_tile(aes(Strain,name,fill=value2), color = "gray70")+
  guides(fill=guide_legend(nrow = 1))+
  scale_fill_manual(name="log10 titer (PFU/mL)", values = my_colors2)+ 
  facet_grid(.~ST_KL, scales="free_x", space="free", labeller = as_labeller(ST_KL_names))+
  theme(legend.box = "horizontal",
    legend.position = "top",
    legend.title=element_text(size=5),
    legend.text=element_text(size=5),
    legend.key.size = unit(3, "mm"),
    legend.margin=margin(b = -5),
    axis.title = element_blank(),
    #axis.text.x =element_text(angle=90, hjust=0, vjust=0.5, size = 4),#
    axis.text.x =element_blank(),
    axis.text.y =element_text(size = 5),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    strip.text = element_text(size=5, margin = unit(c(0,0,0,0), "cm")),
    strip.background = element_rect(color="gray70", fill="white"),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.spacing.x = unit(0.2,"line"),
    plot.margin = unit(c(0,0,0,0), "cm"))+
  scale_x_discrete(position="bottom")+
  scale_y_discrete(limits=rev, position = "left")

plot(plot1)

# include columns for plot 2 and 3

ACI_long_HWM <- ACI_long %>% filter(name == "Highwayman")
ACI_long_HWM <- mutate(ACI_long_HWM, Countries = case_when(name == "Highwayman" ~ "Country"))
ACI_long_HWM <- mutate(ACI_long_HWM, Cities = case_when(name == "Highwayman" ~ "City"))

tbl_cities <- ACI_long_HWM %>% 
  group_by(country, city) %>% 
  summarise(n = n())

ACI_long_HWM$country_city <- paste(ACI_long_HWM$country, "-", ACI_long_HWM$city)

ACI_long_HWM %>% distinct(country_city)


plot2 <- ACI_long_HWM %>% 
  mutate(Strain=factor(Strain, ordered = TRUE,levels = order_Strain) ) %>% 
  mutate(ST_KL=factor(ST_KL, ordered = TRUE,levels=c("ST2 - KL3", "ST636 - KL40", "ST492 - KL104", "ST2 - KL2", "ST1 - KL1", "ST2 - KL9" , "ST2 - KL12" , "ST1 - KL17" ,
                                                     "ST2 - KL7" , "ST2 - KL32" , "ST2 - KL77"))) %>%  
  ggplot()+
  geom_tile(aes(Strain,Cities,fill=city))+ 
  scale_fill_manual(values = c("#30123BFF", "#3D358CFF", "#4457C8FF", "#4777EFFF", "#4195FFFF", "#2EB3F3FF", "#1BD0D5FF", "#1AE4B6FF", "#35F393FF", "#62FC6BFF",
                               "#90FF48FF" ,"#B3F836FF", "#D2E935FF", "#EBD339FF", "#FABA39FF", "#FE9B2DFF", "#F9771EFF", "#EC5410FF", "#DB3A07FF", "#C22402FF",
                                "#A11201FF" ,"#7A0403FF"),
                     name = "")+
  #scale_fill_manual(values = c("forestgreen","chartreuse1","seagreen2","darkolivegreen1","yellowgreen",
                              # "chartreuse3", "olivedrab", "darkseagreen3", "seagreen3",
                              # "royalblue", "royalblue4", "dodgerblue",
                              # "cyan1", "cyan3", "dodgerblue3", "lightblue2",
                              # "yellow1", "gold2", "yellow2",
                              # "purple1", "purple4",
                               #'red1'),
                   # breaks = c("Hungary - szeged","Hungary - budapest", "Hungary - karcag", "Hungary - nyíregyháza", "Hungary - kisvárda",
                              # "Hungary - győr", "Hungary - salgótarján", "Hungary - debrecen", "Hungary - miskolc", 
                               #"Romania - targu_mures","Romania - bacau", "Romania - satu_mare",
                               #"Romania - timisoara","Romania - bucharest","Romania - cluj_napoca","Romania - odorheiu_secuiesc",
                                #"Serbia - belgrade", "Serbia - kragujevac","Serbia - uzice",
                               #"Bosnia and Herzegovina - prijedor","Bosnia and Herzegovina - banja_luka", 
                              # "Montenegro - podgorica"), name = "")+#
  guides(colour=guide_legend("", override.aes=list(fill="azure3", colour="transparent")))+
  ggh4x::facet_nested(.~ST_KL, scales="free_x", space="free")+
  ylab("")+
  xlab(NULL)+
  theme(axis.text.x = element_blank(),
        axis.text.y =element_text(size = 5),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_rect(color="gray70", fill="white"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(color="transparent", fill="transparent"),
        panel.spacing.x = unit(1,"mm"),
        plot.background = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  scale_y_discrete(position = "left")



plot3 <- ACI_long_HWM %>% 
  mutate(Strain=factor(Strain, ordered = TRUE,levels = order_Strain) ) %>% 
  mutate(ST_KL=factor(ST_KL, ordered = TRUE,levels=c("ST2 - KL3", "ST636 - KL40", "ST492 - KL104", "ST2 - KL2", "ST1 - KL1", "ST2 - KL9" , "ST2 - KL12" , "ST1 - KL17" ,
                                                     "ST2 - KL7" , "ST2 - KL32" , "ST2 - KL77"))) %>%  
  ggplot()+
  geom_tile(aes(Strain,Countries,fill=country))+  
  scale_fill_manual(values = c("#C89E6D", "#7ABB9B", "#7590CF","#ffc600", "#7fffdf"),
                    breaks = c("Hungary","Romania","Serbia", "Bosnia and Herzegovina","Montenegro"),
                    name = "")+
  guides(colour=guide_legend("", override.aes=list(fill="azure3", colour="transparent")))+
  ggh4x::facet_nested(.~ST_KL, scales="free_x", space="free")+
  ylab(NULL)+
  xlab(NULL)+
  theme(axis.text.x = element_blank(),
        axis.text.y =element_text(size = 5),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_rect(color="gray70", fill="white"),
        legend.position = "bottom",
        legend.text=element_text(size=5),
        legend.key.size = unit(3, "mm"),
        legend.margin=margin(t=-10),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(1,"mm"),
        plot.background = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  scale_y_discrete(position = "left")
        


#plot(plot3)

#g <- plot1/plot2/plot3 + plot_layout(heights = c(0.5, 0.03, 0.03, 0.1),)

#g <- plot1/plot3 + plot_layout(heights = c(0.5, 0.02),)

g <- plot1/plot_spacer()/plot2/plot_spacer()/plot3 + 
  plot_layout(heights = c(0.5, -0.09, 0.02, -0.09, 0.02),)

plot(g)

ggsave(
  filename = "Fig3A_heatmap_phages_PFU_composed_country_order_noLfw_cities_old_coloring.pdf",
  plot = g,
  width = 180,
  height = 60,
  dpi = 800, units = "mm"
)

saveRDS(g, "Fig3A_heatmap_phages_PFU_composed_country_order_noLfw_cities_old_coloring_pdf.rds")

#g2 <- plot1/plot_spacer()/plot2/plot_spacer()/plot3 + plot_layout(heights = c(0.5, -0.115, 0.02, -0.115, 0.02))

ggsave(
  filename = "Fig3A_heatmap_phages_PFU_composed_country_order_noLfw_cities_old_coloring.png",
  plot = g,
  width = 180,
  height = 60,
  dpi = 800, units = "mm")
  
saveRDS(g, "Fig3A_heatmap_phages_PFU_composed_country_order_noLfw_cities_old_coloring_png.rds")

