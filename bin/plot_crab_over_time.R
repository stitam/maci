rm(list = ls())

library(dplyr)
library(ggplot2)

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  aci_all_path <- args[1]
} else {
  test_dir <- "~/Methods/aci-tests/test-plot_crab_over_time"
  aci_all_path <- paste0(test_dir, "/aci_filtered.rds")
}

aci_all <- readRDS(aci_all_path) %>% dplyr::filter(filtered == TRUE)

# ensure all isolates pass QC
if (any(aci_all$qc_pass == FALSE)) {
  stop("Some isolates do not pass QC. Check.")
}

# if collection_year is below 1945, set to NA
aci_all$collection_year[which(aci_all$collection_year < 1945)] <-  "NA"

# if collection_year is below 2000, pool to <2000
aci_all$collection_year[which(aci_all$collection_year < 2000)] <-  "<2000"

aci <- aci_all %>% 
  filter(collection_year!="NA")  %>%  
  group_by(collection_year, crab) %>% 
  summarise(count=length(assembly))

# ensure first category on plot is <2000
levels <- c(NA, unique(aci$collection_year))
levels[1] <- levels[length(levels)]
levels <- levels[1:(length(levels) -1)]

aci$collection_year <- factor(aci$collection_year, levels = levels)

g <- ggplot(aci, aes(collection_year, count, fill=crab))+
  geom_bar(stat="identity", position=position_dodge(), width = 0.8)+
  xlab("Year of collection")+
  ylab("Number of sequenced isolates")+
  theme_classic()+
  scale_y_continuous(minor_breaks = seq(0 , 2000, 100))+
  scale_fill_manual(
    values=c("gray50","seagreen4"),
    labels=c("Non-CRAB", "CRAB"),
    name=""
  )+
  # guides(fill = "none") + 
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 5),
    axis.title.x = element_text(
      margin = margin(t = 0, r = 0, b = 0, l = 0)
    ),
    axis.title.y = element_text(
      margin = margin(t = 0, r = 0, b = 0, l = 0)
    ),
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    panel.grid.major.y=element_line(color="gray", linetype = 3, linewidth = 0.1),
    axis.line=element_line(color="black", linewidth=0.1),
    axis.ticks = element_line(linewidth=0.1),
    legend.position = c(0.2,0.2),
    legend.key.size = unit(2, "mm"),
    plot.margin = unit(c(0,0,0,0), "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )

ggsave(
    "crab_over_time.pdf", 
    plot = g
)

ggsave(
  "crab_over_time.png", 
  plot = g,
)

saveRDS(g, file = "crab_over_time.rds")