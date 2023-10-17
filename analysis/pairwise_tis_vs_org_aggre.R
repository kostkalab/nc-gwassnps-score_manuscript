library(PRROC)
library(pROC)
library(GenomicRanges)
library(readr)

DHS_SNPSNAP_logreg = readRDS("./data/analysis/bootstrap_DHS_SNPSNAP_logreg.RDS")

bootstrap_org_auc = readRDS("./data/analysis/bootstrap_org_auc.RDS")

DHS_min = DHS_SNPSNAP_logreg[DHS_SNPSNAP_logreg$regularization == "min",]

# get the mean
DHS = tapply(DHS_min$auc, list(DHS_min$PTPE, DHS_min$measure), mean)

bootstrap_org_auc = do.call("rbind", bootstrap_org_auc)
org_auc_mean =  tapply(bootstrap_org_auc$auc, 
                       list(bootstrap_org_auc$PTPE, bootstrap_org_auc$measure, bootstrap_org_auc$variant_score), 
                       mean)

get_p_value = function(org_name, org_auc, measure){
  org_auc = org_auc_mean[,,org_name][,measure]
  tis_auc = DHS[,measure]
  a = cbind(org_auc, tis_auc)
  res = wilcox.test(a[,1], a[,2], paired = T)
  p_val = res$p.value
  score_median = c(median(a[,1]), median(a[,2]))
  names(score_median) = c(org_name, "DHS-Weighted")
  p_df = data.frame(score1 = org_name, score1_median = score_median[1], 
                    score2 = "DHS-Weighted", score2_median = score_median[2], 
                    measure = measure, p_value = p_val, 
                    higher_median = names(which.max(score_median)))
  return(p_df)
}

os = as.character(unique(bootstrap_org_auc$variant_score))
p1 = lapply(os, get_p_value, org_auc_mean, "pr")
p2 = lapply(os, get_p_value, org_auc_mean, "roc")

p1 = do.call("rbind", p1)
p2 = do.call("rbind", p2)
rownames(p1) = NULL
rownames(p2) = NULL
p = rbind(p1, p2)

p %>% dplyr::mutate(dplyr::across(c(2,4), format, digits = 3))  %>% 
  dplyr::mutate(dplyr::across(p_value, format, digits = 4))  %>% 
  readr::write_csv(file="./sup_data/sup_data_pairwise-tis-vs-org-aggregated.csv.gz") 


