library(readr)

geno_SNPSNAP_logreg = readRDS("./data/analysis/bootstrap_geno_SNPSNAP_logreg.RDS")
DHS_SNPSNAP_logreg = readRDS("./data/analysis/bootstrap_DHS_SNPSNAP_logreg.RDS")
Fit_SNPSNAP_logreg = readRDS("./data/analysis/bootstrap_Fit_SNPSNAP_logreg.RDS")
all_SNPSNAP_logreg = rbind(geno_SNPSNAP_logreg, DHS_SNPSNAP_logreg, Fit_SNPSNAP_logreg)

tis_min = all_SNPSNAP_logreg[all_SNPSNAP_logreg$regularization == "min",]
tis_aggre = tapply(tis_min$auc, list(tis_min$PTPE, tis_min$variant_score, tis_min$measure), mean)

tis_base       = all_SNPSNAP_logreg[all_SNPSNAP_logreg$regularization == "baseline",]
tis_aggre_base = tapply(tis_base$auc, list(tis_base$PTPE, tis_base$variant_score, tis_base$measure), mean)

vs = c("DHS", "Fitcons2", "Genoskyline")
com = as.data.frame(combn(vs,2))

get_p_value = function(x, aggre, measure){
  x = as.character(x)
  a = aggre[,,measure][,x]
  res = wilcox.test(a[,1], a[,2], paired = T)
  p_val = res$p.value
  score_median = c(median(a[,x[1]]), median(a[,x[2]]))
  names(score_median) = x
  p_df = data.frame(score1 = x[1], score1_median = score_median[1], 
                    score2 = x[2], score2_median = score_median[2], 
                    measure = measure, p_value = p_val, 
                    higher_median = names(which.max(score_median)))
  rownames(p_df) = NULL
  return(p_df)
}

p_value_pr = lapply(com, get_p_value, tis_aggre, measure = "pr")
p_value_roc = lapply(com, get_p_value, tis_aggre, measure = "roc")

p_value_pr = do.call("rbind", p_value_pr)
p_value_roc = do.call("rbind", p_value_roc)
p_value = rbind(p_value_pr, p_value_roc)
rownames(p_value) = NULL

p_value %>% dplyr::mutate(dplyr::across(c(2,4), format, digits = 4))  %>% 
  dplyr::mutate(dplyr::across(p_value, format, digits = 4))  %>% 
  dplyr::mutate(dplyr::across(c(1,3,7), paste0, "-Weighted")) %>%  # add weighted at the end to let people know it's the logreg
  readr::write_csv(file="./sup_data/sup_data_pairwise-tis-weighted-aggregated.csv.gz") 



get_p_value_paired = function(name, b1 = tis_aggre, b2 = tis_aggre_base, measure = "pr"){
  baseline = b2[,,measure][,name]
  logreg   = b1[,,measure][,name]
  a = cbind(baseline, logreg)
  res = wilcox.test(a[,1], a[,2], paired = T)
  p_val = res$p.value
  score_median = c(median(a[,1]), median(a[,2]))
  names(score_median) = c("Mean", "Weighted")
  p_df = data.frame(score = name, 
                    score1 = "Weighted", score1_median = score_median[1], 
                    score2 = "Mean", score2_median = score_median[2], 
                    measure = measure, p_value = p_val, 
                    higher_median = names(which.max(score_median)))
  rownames(p_df) = NULL
  return(p_df)
}

p_value_pr_paired  = lapply(vs, get_p_value_paired, measure = "pr")
p_value_roc_paired = lapply(vs, get_p_value_paired, measure = "roc")
p_value_pr_paired = do.call("rbind", p_value_pr_paired)
p_value_roc_paired = do.call("rbind", p_value_roc_paired)
p_value_paired = rbind(p_value_pr_paired, p_value_roc_paired)



p_value_paired$score1 = paste(p_value_paired$score, p_value_paired$score1, sep = "-")
p_value_paired$score2 = paste(p_value_paired$score, p_value_paired$score2, sep = "-")
p_value_paired$higher_median = paste(p_value_paired$score, p_value_paired$higher_median, sep = "-")

p_value_paired[,2:8] %>% dplyr::mutate(dplyr::across(c(2,4), format, digits = 4))  %>% 
  dplyr::mutate(dplyr::across(p_value, format, digits = 4))  %>% 
  readr::write_csv(path="./sup_data/sup_data_pairwise-tis-aggregated.csv.gz") 

