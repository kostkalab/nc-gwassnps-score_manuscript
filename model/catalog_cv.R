library(pROC)
library(ggplot2)
library(reshape2)
library(glmnet)
library(pamr)
library(ggbeeswarm)
library(xgboost)
library(PRROC)


source("./model/fun_logreg_cv.R")


# work with the file and all phenotypes
multiPheno = function(file, times){
  dat = readRDS(file)
  col         = colnames(dat)[grep("E[0-9]+", colnames(dat))]
  dat = dat[!is.na(rowSums(dat[,col])),] # remove the rows with all NAs
  dat$phenotype = tolower(gsub(" ","_",dat$phenotype)) # remove space and make lower case
  
  PTPE = unique(dat$phenotype) # go through all phenotypes
  multi_p = lapply(PTPE, rep_procPheno, dat = dat, times = times)
  multi_p = do.call("rbind", multi_p)
  return(multi_p)
}


## read in file
geno_SNPSNAP = "./data/anno_genoskyline_catalog_cutoff100_snpsnap.RDS"
DHS_SNPSNAP  = "./data/anno_DHS_catalog_cutoff100_snpsnap.RDS"
Fit_SNPSNAP  = "./data/anno_Fitcons2_catalog_cutoff100_snpsnap.RDS"



## do logreg
geno_SNPSNAP_logreg = multiPheno(geno_SNPSNAP, times = 30)
DHS_SNPSNAP_logreg  = multiPheno(DHS_SNPSNAP,  times = 30)
Fit_SNPSNAP_logreg  = multiPheno(Fit_SNPSNAP,  times = 30)


geno_SNPSNAP_logreg$variant_score = "Genoskyline"
DHS_SNPSNAP_logreg$variant_score = "DHS"
Fit_SNPSNAP_logreg$variant_score = "Fitcons2"


saveRDS(geno_SNPSNAP_logreg, "./data/analysis/perf_geno_logreg.RDS")
saveRDS(DHS_SNPSNAP_logreg, "/data/analysis/perf_DHS_logreg.RDS")
saveRDS(Fit_SNPSNAP_logreg, "./data/analysis/perf_Fit_logreg.RDS")