library(PRROC)
library(pROC)
library(dplyr)

auc_all = readRDS("./data/analysis/perf_all_divan.RDS")


names(auc_all) = c("measure", "id", "auc", "PTPE", "variant_score")
auc_all[auc_all$variant_score == "min","variant_score"] = "DHS_logreg"


auc_split = split(auc_all, f = auc_all$PTPE)

get_p_value = function(x, measure){
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


p_value_pr = lapply(auc_split, get_p_value, measure = "pr")
p_value_roc = lapply(auc_split, get_p_value, measure = "roc")
p_value = rbind(do.call("rbind", p_value_pr), do.call("rbind", p_value_roc))

# pick the top and second one
top_p = lapply(p_value_pr, function(x){
  higher = as.character(unique(x$higher_median))
  p = x[x$score1 %in% higher & x$score2 %in% higher,]
  return(p)
})

top_p = do.call("rbind", top_p)
saveRDS(top_p, "./data/analysis/pairwise_DIVAN_top.RDS")


# suppl.15
p_value  %>% dplyr::mutate(dplyr::across(c(2,4), format, digits = 3))  %>% 
  dplyr::mutate(dplyr::across(p_value, format, digits = 4))  %>% 
  readr::write_csv(file ="./sup_data/sup_data_pairwise-divan-individual.csv.gz") 






