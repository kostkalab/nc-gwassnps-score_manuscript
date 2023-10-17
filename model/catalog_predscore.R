library(pROC)
library(ggplot2)
library(reshape2)
library(glmnet)
library(pamr)
library(ggbeeswarm)
library(PRROC)
library(readr)


source("./model/fun_logreg_full_predscore.R")

multi_predscore = function(file){
  dat = readRDS(file)
  col         = colnames(dat)[grep("E[0-9]+", colnames(dat))]
  dat = dat[!is.na(rowSums(dat[,col])),]
  dat$phenotype = tolower(gsub(" ","_",dat$phenotype))
  
  PTPE = unique(dat$phenotype)
  
  ## get the prediction score and the model
  pred_all = lapply(PTPE, get_prediction_score, dat = dat)
  # get prediction scores
  pred_score_all = lapply(pred_all, function(x) x[[2]])
  pred_score_all = do.call("rbind", pred_score_all)
  # get the model
  model_all = lapply(pred_all, function(x) x[[1]])
  names(model_all) = PTPE
  
  return(list(pred_score_all, model_all))
}


## read in file
geno_SNPSNAP = "./data/anno_genoskyline_catalog_cutoff100_snpsnap.RDS"
DHS_SNPSNAP  = "./data/anno_DHS_catalog_cutoff100_snpsnap.RDS"
Fit_SNPSNAP  = "./data/anno_Fitcons2_catalog_cutoff100_snpsnap.RDS"

## get prediction score
geno_pred = multi_predscore(geno_SNPSNAP)
DHS_pred  = multi_predscore(DHS_SNPSNAP)
Fit_pred  = multi_predscore(Fit_SNPSNAP)


# prediction scores
predic_score = cbind(geno_pred[[1]], DHS_pred[[1]][,4], Fit_pred[[1]][,4])
colnames(predic_score)[4:6] = c("Genoskyline_Weighted", "DHS_Weighted", "Fitcon2_Weighted")
colnames(predic_score)[1] = "SNV_ID"
readr::write_csv(predic_score, path="./sup_data/sup_data_prediction-scores-dhs-weighted.csv.gz") 

# save model
saveRDS(geno_pred[[2]], "./data/analysis/model_111diseases_Genoskyline.RDS")
saveRDS(DHS_pred[[2]], "./data/analysis/model_111diseases_DHS.RDS")
saveRDS(Fit_pred[[2]], "./data/analysis/model_111diseases_Fitcons2.RDS")
