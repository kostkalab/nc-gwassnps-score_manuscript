library(PRROC)
library(pROC)
library(GenomicRanges)
library(dplyr)
library(readr)

source("./analysis/fun_bootstrap_pairwise_org.R")

catalog_org         <- readRDS("./data/anno_org_catalog_cutoff100_snpsnap.RDS")

disease = catalog_org[catalog_org$type == 1]
sc = c("CADD", "eigen", "GenoCanyon", "GWAVA", "LINSIGHT")
t = table(catalog_org[catalog_org$type == 1]$phenotype)
phenotypes = sort(names(t[t>=100]))


bootstrap_org_auc = lapply(phenotypes, bootstrap_auc, org = catalog_org, sc = sc) # why cannot use GenoCanyon? 
saveRDS(bootstrap_org_auc, "./data/analysis/bootstrap_org_auc.RDS")
