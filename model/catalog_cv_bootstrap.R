library(GenomicRanges)

source("./model/fun_logreg_cv.R")

# use two functions from this 
# mk_model
# procPhenotype


library(sampler)
library(dplyr)
library(PRROC)
library(pROC)

bootstrap_samples = readRDS("./data/analysis/bootstrap_samples_30_0.9.RDS")
names(bootstrap_samples) = tolower(gsub(" ","_",names(bootstrap_samples)))

bootstrap_glmnet = function(PTPE, dat){
  set.seed(0525) 
  dat_PTPE = dat[dat$phenotype == PTPE,]
  bootstrap_PTPE = data.frame(bootstrap_samples[[PTPE]])
  colnames(bootstrap_PTPE) = 1:30
  res_cd = lapply(bootstrap_PTPE[,1:30], function(x){
    dat_partial = dat_PTPE[dat_PTPE$snpID %in% x, ]
    res = procPhenotype(PTPE, dat_partial)
    return(res)
  })
  df = melt(res_cd); colnames(df) = c("regularization","measure","auc","id")
  df$PTPE = PTPE
  return(df)
}


bootstrapPheno = function(file){
  dat = readRDS(file)
  dat = as.data.frame(dat)
  col = colnames(dat)[grep("E[0-9]+", colnames(dat))]
  dat = dat[!is.na(rowSums(dat[,col])),]
  dat$phenotype = tolower(gsub(" ","_",dat$phenotype))

  PTPE = unique(dat$phenotype)
  multi_p = lapply(PTPE, bootstrap_glmnet, dat = dat)
  multi_p = do.call("rbind", multi_p)
  return(multi_p)
}


## read in file
geno_SNPSNAP = "./data/anno_genoskyline_catalog_cutoff100_snpsnap.RDS"
DHS_SNPSNAP  = "./data/anno_DHS_catalog_cutoff100_snpsnap.RDS"
Fit_SNPSNAP  = "./data/anno_Fitcons2_catalog_cutoff100_snpsnap.RDS"


## do bootstrap
geno_SNPSNAP_logreg = bootstrapPheno(geno_SNPSNAP)
DHS_SNPSNAP_logreg  = bootstrapPheno(DHS_SNPSNAP)
Fit_SNPSNAP_logreg  = bootstrapPheno(Fit_SNPSNAP)



geno_SNPSNAP_logreg$variant_score = "Genoskyline"
DHS_SNPSNAP_logreg$variant_score = "DHS"
Fit_SNPSNAP_logreg$variant_score = "Fitcons2"


saveRDS(geno_SNPSNAP_logreg, "./data/analysis/bootstrap_geno_SNPSNAP_logreg.RDS")
saveRDS(DHS_SNPSNAP_logreg,  "./data/analysis/bootstrap_DHS_SNPSNAP_logreg.RDS")
saveRDS(Fit_SNPSNAP_logreg, "./data/analysis/bootstrap_Fit_SNPSNAP_logreg.RDS")