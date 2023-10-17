library(dplyr)
library(readr)

bootstrap_org_auc = readRDS("./data/analysis/bootstrap_org_auc.RDS")

get_p_value_org_disease = function(x, measure){
  org_metric = x[x$measure == measure, ]
  vs = as.character(unique(org_metric$variant_score))
  com = as.data.frame(combn(vs,2))
  p = lapply(com, function(x){
    x = as.character(x)
    print(x)
    a1 = org_metric[org_metric$variant_score == x[1], c("id", "auc")]
    colnames(a1)[2] = x[1]
    a2 = org_metric[org_metric$variant_score == x[2], c("id", "auc")]
    colnames(a2)[2] = x[2]
    a = full_join(a1, a2, by = "id")
    res = wilcox.test(a[,2], a[,3], paired = T)
    p_val = res$p.value
    score_median = c(median(a[,2]), median(a[,3]))
    names(score_median) = x
    p_df = data.frame(score1 = x[1], score1_median = score_median[1], 
                      score2 = x[2], score2_median = score_median[2], 
                      disease = org_metric$PTPE[1], 
                      measure = measure, p_value = p_val, 
                      higher_median = names(which.max(score_median)))
    return(p_df)
  })
  p_value = do.call("rbind", p)
  rownames(p_value) = NULL
  return(p_value)
}


p_value_disease_roc = lapply(bootstrap_org_auc, get_p_value_org_disease, measure = "roc")
p_value_disease_pr = lapply(bootstrap_org_auc, get_p_value_org_disease, measure = "pr")

p_value_disease_roc = do.call("rbind", p_value_disease_roc)
p_value_disease_pr = do.call("rbind", p_value_disease_pr)

p_value_disease = rbind(p_value_disease_roc, p_value_disease_pr)
#p_value_disease %>% dplyr::mutate(dplyr::across(where(is.numeric), format, digits = 4)) %>% 
#   readr::write_csv(path="~/projects/variant-scores-2/data/GWAS_catalog_project_3.0/genetic/suppl4.csv") 
p_value_disease %>% dplyr::mutate(dplyr::across(c(2,4), format, digits = 3))  %>% 
  dplyr::mutate(dplyr::across(p_value, format, digits = 4))  %>% 
   readr::write_csv(file="./sup_data/sup_data_pairwise-org-individual.csv.gz") 


