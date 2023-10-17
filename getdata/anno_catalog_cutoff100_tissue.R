library(GenomicRanges)
library(seqminer)
library(rtracklayer)
library(ggplot2)
source("./getdata/fun_read_variantscore_tissue.R")

matched_data2 = readRDS("./data/catalog_cutoff100_snpsnap.RDS")


genoskyline_SNPSNAP2 = read_genoskyline(matched_data2)
DHS_SNPSNAP2         = read_DHS(matched_data2)
Fitcons2_SNPSNAP2    = read_Fitcons2(matched_data2)


saveRDS(genoskyline_SNPSNAP2, "./data/anno_genoskyline_catalog_cutoff100_snpsnap.RDS")
saveRDS(DHS_SNPSNAP2,         "./data/anno_DHS_catalog_cutoff100_snpsnap.RDS")
saveRDS(Fitcons2_SNPSNAP2,    "./data/anno_Fitcons2_catalog_cutoff100_snpsnap.RDS")
