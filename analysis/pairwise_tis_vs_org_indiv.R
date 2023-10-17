library(PRROC)
library(pROC)
library(GenomicRanges)


DHS_SNPSNAP_logreg = readRDS("./data/analysis/bootstrap_DHS_SNPSNAP_logreg.RDS")

bootstrap_org_auc = readRDS("./data/analysis/bootstrap_org_auc.RDS")

# tis pick DHS
# organism pick genocanyon


DHS_min = DHS_SNPSNAP_logreg[DHS_SNPSNAP_logreg$regularization == "min",]


get_p_value = function(x, measure){
  phenotype = tolower(gsub(" ","_",x$PTPE[1]))
  
  or = as.character(unique(x$variant_score))
  p = lapply(or, function(or_s){
    org = x[x$variant_score == or_s & x$measure == measure, c("auc", "id")]
    colnames(org)[1] = or_s
    tis = DHS_min[DHS_min$PTPE == phenotype & DHS_min$measure == measure, c("auc", "id")]
    tis$id = as.integer(tis$id)
    colnames(tis)[1] = "DHS"
    com = full_join(org, tis, by = "id")
    res = wilcox.test(com[,or_s], com[,"DHS"], paired = T)
    p_val = res$p.value
    score_median = c(median(com[,or_s]), median(com[,"DHS"]))
    names(score_median) = c(or_s, "DHS-Weighted")
    p_df = data.frame(score1 = or_s, score1_median = score_median[1], 
                      score2 = "DHS-Weighted", score2_median = score_median[2], 
                      disease = phenotype, 
                      measure = measure, p_value = p_val, 
                      higher_median = names(which.max(score_median)))
    return(p_df)
  })
  p_value = do.call("rbind", p)
  rownames(p_value) = NULL
  return(p_value)
}

p_value_pr  = lapply(bootstrap_org_auc, get_p_value, measure = "pr")
p_value_roc = lapply(bootstrap_org_auc, get_p_value, measure = "roc")

p_value_pr2  = do.call("rbind", p_value_pr)
p_value_roc2 = do.call("rbind", p_value_roc)
p_value = rbind(p_value_pr2, p_value_roc2)

p_value %>% dplyr::mutate(dplyr::across(c(2,4), format, digits = 3))  %>% 
  dplyr::mutate(dplyr::across(p_value, format, digits = 4))  %>% 
  readr::write_csv(path="./sup_data/sup_data_pairwise-tis-vs-org-individual.csv.gz") 

