require(lpSolve)
library(dplyr)

source("./getdata/fun_chrom_heldout_annotation.R")

snpsnap_100 = readRDS("./data/catalog_cutoff100_snpsnap.RDS")
phenotypes = unique(snpsnap_100$phenotype)

snpsnap_100_chrom = lapply(phenotypes, get_chrom_heldout, snpsnap_100)
snpsnap_100_chrom = do.call("c", snpsnap_100_chrom)

saveRDS(snpsnap_100_chrom, "./data/catalog_cutoff100_snpsnap_chrom_heldout.RDS")
