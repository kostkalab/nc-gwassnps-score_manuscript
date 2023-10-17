library(pROC)
library(ggplot2)
library(reshape2)
library(glmnet)
library(pamr)
library(ggbeeswarm)
library(xgboost)
library(PRROC)
library(dplyr)


source("./model/fun_logreg_cv.R")


## read in file
geno_SNPSNAP = "./data/anno_genoskyline_catalog_cutoff100_snpsnap.RDS"
DHS_SNPSNAP  = "./data/anno_DHS_catalog_cutoff100_snpsnap.RDS"
Fit_SNPSNAP  = "./data/anno_Fitcons2_catalog_cutoff100_snpsnap.RDS"

ldfree_marked <- readRDS("./data/catalog_cutoff100_snpsnap_ld_marked.RDS")
ldfree_marked$phenotype =  tolower(gsub(" ","_",ldfree_marked$phenotype))
ldfree_marked$id = paste(ldfree_marked$snpID, ldfree_marked$phenotype, sep = "_") # generate a id column, so each snp id is unique among all. use it to match


multiPheno_ld = function(file, times){
  dat = readRDS(file)
  
  col         = colnames(dat)[grep("E[0-9]+", colnames(dat))]
  dat = dat[!is.na(rowSums(dat[,col])),] # remove rows with all NAs
  dat$phenotype = tolower(gsub(" ","_",dat$phenotype)) # remove space and all to lower case
  
  dat$id = paste(dat$snpID, dat$phenotype, sep = "_")
  
  a = left_join(dat, as.data.frame(mcols(ldfree_marked))[,c("id", "ld_block_picked")], by = "id") # get the ld_block_picked column for dat
  dat = a[a$ld_block_picked == 1, !colnames(a) %in% "id"] # only keep picked snps and delete the id column which will not be used later
  PTPE = unique(dat$phenotype)
  
  ## get the auc of cross validation
  multi_p = lapply(PTPE, rep_procPheno, dat = dat, times = times)
  multi_p = do.call("rbind", multi_p)
  
  return(multi_p)
}


## do logreg
geno_SNPSNAP_logreg = multiPheno_ld(geno_SNPSNAP, times = 30)
DHS_SNPSNAP_logreg  = multiPheno_ld(DHS_SNPSNAP,  times = 30)
Fit_SNPSNAP_logreg  = multiPheno_ld(Fit_SNPSNAP,  times = 30)


geno_SNPSNAP_logreg$variant_score = "Genoskyline"
DHS_SNPSNAP_logreg$variant_score = "DHS"
Fit_SNPSNAP_logreg$variant_score = "Fitcons2"


saveRDS(geno_SNPSNAP_logreg, "./data/analysis/perf_geno_logreg_ldfree.RDS")
saveRDS(DHS_SNPSNAP_logreg, "./data/analysis/perf_DHS_logreg_ldfree.RDS")
saveRDS(Fit_SNPSNAP_logreg, "./data/analysis/perf_Fit_logreg_ldfree.RDS")