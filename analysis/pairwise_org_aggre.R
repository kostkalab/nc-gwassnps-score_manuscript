library(PRROC)
library(pROC)
library(GenomicRanges)
library(tidyverse)
catalog_org_auc <- readRDS("./data/analysis/perf_org.RDS")
catalog_org_auc$PTPE = as.character(catalog_org_auc$PTPE)
catalog_org_auc = catalog_org_auc[catalog_org_auc$match == "SNPsnap",]

catalog_org_auc = catalog_org_auc[,1:4]

get_p_value_org = function(x, measure){
  org_metric = x[x$measure == measure, ]
  vs = as.character(unique(org_metric$variant_score))
  com = as.data.frame(combn(vs,2))
  p = lapply(com, function(x){
    x = as.character(x)
    print(x)
    a1 = org_metric[org_metric$variant_score == x[1], c("PTPE", "auc")]
    colnames(a1)[2] = x[1]
    a2 = org_metric[org_metric$variant_score == x[2], c("PTPE", "auc")]
    colnames(a2)[2] = x[2]
    a = full_join(a1, a2, by = "PTPE")
    res = wilcox.test(a[,2], a[,3], paired = T)
    p_val = res$p.value
    score_median = c(median(a[,2]), median(a[,3]))
    names(score_median) = x
    p_df = data.frame(score1 = x[1], score1_median = score_median[1], 
                      score2 = x[2], score2_median = score_median[2], 
                      measure = measure, p_value = p_val, 
                      higher_median = names(which.max(score_median)))
    return(p_df)
  })
  p_value = do.call("rbind", p)
  rownames(p_value) = NULL
  return(p_value)
}

p_value_pr = get_p_value_org(catalog_org_auc, "pr")
p_value_roc  = get_p_value_org(catalog_org_auc, "roc")

p_value = rbind(p_value_roc, p_value_pr)

p_value %>% dplyr::mutate(dplyr::across(c(2,4), format, digits = 4))  %>% 
  dplyr::mutate(dplyr::across(p_value, format, digits = 4))  %>% 
  readr::write_csv(file = "./sup_data/sup_data_pairwise-org-aggregated.csv.gz") 

