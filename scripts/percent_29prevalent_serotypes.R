library(tidyverse)

tbl_r <- read_tsv("top_serotypes_region23_collapse_geodate.tsv")

tbl_r %>% distinct(serotype)# 223
tbl_r %>% distinct(region23)# 9 

t <- tbl_r %>% group_by(serotype) %>%  
  summarise(Nr = sum(count)) %>%  
  #arrange(desc(Nr))#
  summarise(S = sum(Nr)) # 3150

# the 29 prevalent serotypes
tbl_prev <- tbl_r %>% filter (serotype == "ST1_KL1" | 
                                serotype == "ST1_KL18" |
                                serotype == "ST15_KL22" |
                                serotype == "ST16_KL24"|
                                serotype =="ST2_KL10"|
                                serotype =="ST2_KL12"|
                                serotype =="ST2_KL13"|
                                serotype =="ST2_KL152"|
                                serotype =="ST2_KL160"|
                                serotype =="ST2_KL2"|
                                serotype =="ST2_KL210"|
                                serotype =="ST2_KL22"|
                                serotype =="ST2_KL234"|
                                serotype =="ST2_KL3"|
                                serotype =="ST2_KL4"|
                                serotype =="ST2_KL6"|
                                serotype =="ST2_KL7"|
                                serotype =="ST2_KL77"|
                                serotype =="ST2_KL81"|
                                serotype =="ST2_KL9"|
                                serotype =="ST25_KL139"|
                                serotype =="ST25_KL22"|
                                serotype =="ST3_KL17"|
                                serotype =="ST492_KL104"|
                                serotype =="ST499_KL18"|
                                serotype =="ST636_KL40"|
                                serotype =="ST78_KL3"|
                                serotype =="ST79_KL49"|
                                serotype =="ST79_KL9" ) %>% 
  group_by(serotype) %>%  
  summarise(Nr = sum(count)) %>% 
  summarise(S = sum(Nr)) # 2212

# percent of the prevalent serotypes

p = tbl_prev$S/t$S*100
print(p) # 70,2 %
