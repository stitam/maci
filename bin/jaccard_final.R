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
library(dplyr)
library(reshape2)

# recent dataset (20230810)

###################
###################
# we consider only 3 world regions (#east-sout_europe, western_europe, northern_america)
###################
###################

aci <- readRDS(args$file) # contains all aci strains with metadata
aci <- aci %>% rename(Assembly=assembly)

aci_all <- read.csv(args$sensitivity, sep = "\t") # contains the strains for we have phage sensitivity

#Left  join in R:  merge() function takes df1 and df2 as argument along with all.x=TRUE  
#there by returns all rows from the left table, and any rows with matching keys from the right table.

aci2 <- left_join(aci_all, aci, by ="Assembly")

aci2 %>% distinct(region23) #eastern_europe, NA, southern_europe, western_europe, northern_america

#b <- aci2 %>% filter (is.na (region23)) #101, Aci 198, 212, 215, 298 serotype should be affirmed by sequencing
#c <- aci2 %>% filter (serotype == "NA NA") #97

# introducing new columns

aci2 <- aci2 %>% 
  mutate(serotype = paste(MLST, KL)) %>% 
  group_by(serotype) %>% mutate(id = row_number()) %>% ungroup() %>%  # labels the strains with the same MLST_KL based on the row number
  mutate(serotype2 = paste(serotype, id))

# filtering out NAs

aci2 <- aci2 %>% filter (serotype != "NA NA")
aci2 <- aci2 %>% filter (!is.na (region23))
aci2 <- aci2 %>% filter (if_all( c(Highwayman,Silvergun,Fanak,PhT2.v2,Porter,Dino,Fishpie,Tama,Margaret,ABW132,ABW311,Navy4.v2,Rocket,Konradin.v2,KissB), function(x) !is.na(x)))  

#aci2 %>% distinct(region23)  

# there are 35 strains where none of the phages can infect. From these 17 with different MLST KL.
# I leave those rows in the table where for the respective KL we have phage infection 

no_inf <- aci2 %>% filter (Highwayman == "NO" & Silvergun == "NO" & Fanak == "NO" & PhT2.v2 == "NO" & Porter == "NO" & Dino == "NO" & Fishpie == "NO" & Tama == "NO" & Margaret == "NO" & ABW132 == "NO" &  
                             ABW311 == "NO" & Navy4.v2 == "NO" & Rocket == "NO" & Konradin.v2 == "NO" & KissB == "NO") #%>% group_by(KL) %>% summarise(n = n())

inf <- aci2 %>% filter (Highwayman != "NO" | Silvergun != "NO" | Fanak != "NO" | PhT2.v2 != "NO" | Porter != "NO" | Dino != "NO" | Fishpie != "NO" | Tama != "NO" | Margaret != "NO" | ABW132 != "NO" |  
                          ABW311 != "NO" | Navy4.v2 != "NO" | Rocket != "NO" | Konradin.v2 != "NO" | KissB != "NO") #%>% group_by(KL) %>% summarise(n = n()) #15 KL types

#no_inf_to_keep <- semi_join(no_inf, inf, by ='KL') %>% group_by(KL) %>% summarise(n = n()) # 24 strains, 9 KL types

no_inf_not_to_keep <- anti_join(no_inf, inf, by ='KL') #%>% group_by(KL) %>% summarise(n = n()) # 11 strains, 8 KL types

aci2 <- anti_join(aci2, no_inf_not_to_keep, by ='Assembly')


phage_names <- names(aci2)[8:22]

aci3 <- aci2 %>% select(Assembly, all_of (phage_names), region23, country, serotype, serotype2)

aci3 %>% distinct(region23)

aci3 <- mutate(aci3, regionB3 = case_when(region23 == "eastern_europe" | region23 == "southern_europe" ~ "east_south_europe",
                                          .default = region23))

aci3_matrix <- aci3 %>% select(all_of (phage_names)) %>% as.matrix()
row.names(aci3_matrix) <- aci3$serotype2

aci3_matrix<-toupper(aci3_matrix)
aci3_matrix[aci3_matrix=="1" | aci3_matrix=="NO"]<-"0"
aci3_matrix[aci3_matrix=="YES" ]<-"1"
aci3_matrix[ (as.numeric(aci3_matrix) >5 & as.numeric(aci3_matrix)<12 )]<-"1"
aci3_matrix<-matrix(as.numeric(aci3_matrix) , nrow=nrow(aci3_matrix), dimnames = dimnames(aci3_matrix) ) 


#######################################
############# jaccard index ###########
#######################################


# I kept the same variable names as for euclidean distance script to not rewrite the script

euc_distance_mx<-as.matrix(vegan::vegdist(aci3_matrix, method="jaccard", na.rm=T))


tmp <- aci3 %>%  select(serotype2, regionB3) %>%  distinct() %>%  filter(!is.na(regionB3))


tbl2<-reshape2::melt(euc_distance_mx) %>%
  as_tibble() %>% 
  rename(strain1=Var1, strain2=Var2, euc=value) %>% 
  mutate(
    serotype1=gsub("([^ ]+) ([^ ]+) ([^ ]+)", "\\1 \\2", strain1),
    serotype2=gsub("([^ ]+) ([^ ]+) ([^ ]+)", "\\1 \\2", strain2),
    
    st1=gsub("([^ ]+) ([^ ]+) ([^ ]+)", "\\1", strain1),
    st2=gsub("([^ ]+) ([^ ]+) ([^ ]+)", "\\1", strain2),
    
    kl1=gsub("([^ ]+) ([^ ]+) ([^ ]+)", "\\2", strain1),
    kl2=gsub("([^ ]+) ([^ ]+) ([^ ]+)", "\\2", strain2)
    
  ) %>% 
  inner_join(tmp %>% rename(region_1=regionB3),by = c("strain1"="serotype2") ) %>% 
  inner_join(tmp %>% rename(region_2=regionB3),by = c("strain2"="serotype2") )  %>% 
  filter(strain1!=strain2) %>% 
  group_by(serotype1, serotype2, kl1, kl2, st1, st2, region_1, region_2)  %>%  summarise(
    mean_euc=mean(euc) ,
    n=n()/2 ,
    .groups = "drop")


similarities<-tbl2 %>%
  mutate(
    type=case_when(
      serotype1==serotype2 &  region_1 == region_2 ~ "st_kl_within_same region",
      serotype1==serotype2 &  region_1 != region_2 ~ "st_kl_within_different region",
      kl1==kl2 & st1!=st2 ~"kl_within_st_btw",
      st1 == st2 & kl1 != kl2 ~"st_within_kl_btw",
      st1 != st2 & kl1 != kl2 ~"st_kl_btw",
      .default = "ERROR"
    ), )


tbl_real_means<- similarities %>% group_by(type) %>% summarise(dist=mean(mean_euc,na.rm = TRUE))
tbl_real_mean_diffs<-expand_grid(type1=tbl_real_means$type,type2=tbl_real_means$type) %>% 
  filter(type1>type2) %>% 
  left_join(tbl_real_means %>% rename(dist1=dist,type1=type),by = join_by(type1)) %>% 
  left_join(tbl_real_means %>% rename(dist2=dist,type2=type),by = join_by(type2)) %>%
  mutate(difference_of_distances=dist1-dist2)

similarities$type <- factor(
  similarities$type, 
  levels=c(
    "st_kl_within_same region", 
    "st_kl_within_different region",
    "kl_within_st_btw",
    "st_within_kl_btw",
    "st_kl_btw"
  ),
  ordered = TRUE
)




####################################
######## permutation test ##########
####################################

tmp <- aci3 %>%  select(serotype2, regionB3) %>%  distinct() %>%  filter(!is.na(regionB3))

tbl_bootrstrap_result<-tibble()

# set random seed
set.seed(0)

for(iB in 1:5000){
  euc_distance_mx<-as.matrix(vegan::vegdist(aci3_matrix, method="jaccard", na.rm=T))
  
  
  mxB<-euc_distance_mx
  random_order_names<-sample(rownames(mxB)) #randomizálás
  rownames(mxB)<-random_order_names
  colnames(mxB)<-random_order_names
  
  
  
  tbl2_B<-reshape2::melt(mxB) %>%
    as_tibble() %>% 
    rename(strain1=Var1, strain2=Var2, euc=value) %>% 
    mutate(
      serotype1=gsub("([^ ]+) ([^ ]+) ([^ ]+)", "\\1 \\2", strain1),
      serotype2=gsub("([^ ]+) ([^ ]+) ([^ ]+)", "\\1 \\2", strain2),
      
      st1=gsub("([^ ]+) ([^ ]+) ([^ ]+)", "\\1", strain1),
      st2=gsub("([^ ]+) ([^ ]+) ([^ ]+)", "\\1", strain2),
      
      kl1=gsub("([^ ]+) ([^ ]+) ([^ ]+)", "\\2", strain1),
      kl2=gsub("([^ ]+) ([^ ]+) ([^ ]+)", "\\2", strain2)
      
    ) %>% 
    inner_join(tmp %>% rename(region_1=regionB3),by = c("strain1"="serotype2") ) %>% 
    inner_join(tmp %>% rename(region_2=regionB3),by = c("strain2"="serotype2") )  %>% 
    filter(strain1!=strain2) %>% 
    group_by(serotype1, serotype2, kl1, kl2, st1, st2, region_1, region_2)  %>%  summarise(
      mean_euc=mean(euc) ,
      n=n()/2 ,
      .groups = "drop")
  
  
  similarities_B<-tbl2_B %>%
    mutate(
      type=case_when(
        serotype1==serotype2 &  region_1 == region_2 ~ "st_kl_within_same region",
        serotype1==serotype2 &  region_1 != region_2 ~ "st_kl_within_different region",
        kl1==kl2 & st1!=st2 ~"kl_within_st_btw",
        st1 == st2 & kl1 != kl2 ~"st_within_kl_btw",
        st1 != st2 & kl1 != kl2 ~"st_kl_btw",
        .default = "ERROR"
      ), )
  
  tbl_B2<- similarities_B %>% group_by(type) %>% summarise(dist=mean(mean_euc,na.rm = TRUE))
  tbl_B3<-expand_grid(type1=tbl_B2$type,type2=tbl_B2$type) %>% filter(type1>type2) %>% 
    left_join(tbl_B2 %>% rename(dist1=dist,type1=type),by = join_by(type1)) %>% 
    left_join(tbl_B2 %>% rename(dist2=dist,type2=type),by = join_by(type2)) %>%
    mutate(difference_of_distances=dist1-dist2, i=iB)
  
  tbl_bootrstrap_result<-bind_rows(tbl_bootrstrap_result,tbl_B3)
}

#### histogram with the bootstraped data ####

tbl_bootrstrap_result %>%                                    
  mutate(type_pair=paste(type1,type2)) %>%
  ggplot(mapping =aes(x=difference_of_distances))+
  geom_histogram()+
  geom_vline(data = tbl_real_mean_diffs %>% mutate(type_pair=paste(type1,type2)), mapping = aes(xintercept=difference_of_distances))+
  facet_wrap(.~type_pair)

ggsave("histogram_3regions_bootstrap.pdf", width=10 , height=5)

# making the table with the p values

p_val_tbl <- 
  tbl_bootrstrap_result %>%
  rename(difference_of_distances_rnd=difference_of_distances) %>%
  left_join(tbl_real_mean_diffs, by=c("type1","type2")) %>%
  group_by(type1,type2) %>% 
  summarise(
    p_left=mean(difference_of_distances_rnd<difference_of_distances)
  ) %>% 
  mutate(p=1-2*abs(0.5-p_left)) %>% 
  mutate(signif_label=case_when(p<0.001~"***", p<0.01~"**", p<0.05~"*", .default = "NS")) 

p_val_tbl$p <- format(p_val_tbl$p, scientific = TRUE)

p_val_tbl$shortlisted <- c(
  FALSE,
  TRUE,
  FALSE,
  TRUE,
  FALSE,
  TRUE,
  TRUE,
  TRUE,
  TRUE,
  TRUE
)

p_val_tbl$y <- NA
p_val_tbl$y[which(p_val_tbl$shortlisted)]=c(
  1.55, 1.25, 1.1, 1.25, 1.55, 1.1, 1.4
)

# boxplot with the p-values got from bootstrapping

g <- ggplot(similarities, aes(type, mean_euc))+
  geom_boxplot(outlier.color = "transparent", lwd = 0.1)+
  geom_point(position = position_jitter(0.1, seed = 0), size=0.3, alpha = 0.5)+
  scale_x_discrete(
    name="",
    labels=c(
      "same MLST\nsame CPS\nsame region",
      "same MLST\nsame CPS\ndifferent regions",
      "different MLST\nsame CPS",
      "same MLST\ndifferent CPS",
      "different MLST\ndifferent CPS"
    ),
  )+
  ylim(0, max(p_val_tbl$y, na.rm = TRUE) + 0.1) +
  ylab("Dissimilarity of phage susceptibility profiles \n(Jaccard distance)")+
 
  ggsignif::geom_signif(
    data = p_val_tbl %>% filter(shortlisted==TRUE),
    aes(xmin = type1, 
        xmax = type2, 
        annotations = signif_label, 
        y_position = y),
    textsize = 2,
    tip_length = 0.015,
    manual = TRUE,
    vjust = 0.1,
    size=0.1, map_signif_level=TRUE, step_increase = 0.1, inherit.aes = FALSE)+
  
  theme_classic()+
  theme(axis.text.x  = element_text(size=5, angle = 45, hjust = 1),
        axis.title = element_text(size=5),
        axis.line = element_line(size=0.1))

print(g)

g$pair_counts <- table(similarities$type)
g$p_values <- p_val_tbl

ggsave(
  filename = "Fig3B_sensitivity_dissimilarity.pdf",
  plot = g,
  width = 6,
  height = 6
)
ggsave(
  filename = "Fig3B_sensitivity_dissimilarity.png",
  plot = g,
  width = 6,
  height = 6
)
saveRDS(g, "Fig3B_sensitivity_dissimilarity.rds")

