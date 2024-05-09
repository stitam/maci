downsample_poppunk <- function(aci, pp) {
  
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
  
  return(aci_ds)
  
}