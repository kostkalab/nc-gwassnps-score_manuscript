library(PRROC)
library(pROC)
library(GenomicRanges)
library(dplyr)
library(readr)

# tissue-weighted vs tissue-weighted goes first
# tissue-weighted vs tissue-mean goes next


geno_SNPSNAP_logreg = readRDS("./data/analysis/bootstrap_geno_SNPSNAP_logreg.RDS")
DHS_SNPSNAP_logreg = readRDS("./data/analysis/bootstrap_DHS_SNPSNAP_logreg.RDS")
Fit_SNPSNAP_logreg = readRDS("./data/analysis/bootstrap_Fit_SNPSNAP_logreg.RDS")

all_SNPSNAP_logreg = rbind(geno_SNPSNAP_logreg, DHS_SNPSNAP_logreg, Fit_SNPSNAP_logreg)
all_SNPSNAP_logreg = split(all_SNPSNAP_logreg, f = all_SNPSNAP_logreg$PTPE)

get_p_value_tis_disease = function(x, measure, regularization){
  org_metric = x[x$measure == measure & x$regularization == regularization, ]
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


# baseline vs baseline
baseline_p_value_pr = lapply(all_SNPSNAP_logreg, get_p_value_tis_disease, measure = "pr", regularization = "baseline")
baseline_p_value_roc = lapply(all_SNPSNAP_logreg, get_p_value_tis_disease, measure = "roc", regularization = "baseline")

# logreg vs logreg
logreg_p_value_pr = lapply(all_SNPSNAP_logreg, get_p_value_tis_disease, measure = "pr", regularization = "min")
logreg_p_value_roc = lapply(all_SNPSNAP_logreg, get_p_value_tis_disease, measure = "roc", regularization = "min")


get_stats = function(x){
  names(x) = NULL
  all = do.call("rbind", x)
  all_p = all[all$p_value < 0.05,]
  t = table(all_p$higher_median)
  t1 = as.integer(t)
  names(t1) = names(t)
  t1 = data.frame(t(t1))
  return(t1)
}

s1 = get_stats(baseline_p_value_pr)
s2 = get_stats(baseline_p_value_roc)
s3 = get_stats(logreg_p_value_pr)
s4 = get_stats(logreg_p_value_roc)

s = rbind(s1,s2,s3,s4)
rownames(s) = c("baseline_pr", "baseline_roc", "logreg_pr", "logreg_roc")
s$insig = 333-s$Genoskyline-s$DHS-s$Fitcons2
s



logreg_p_value_pr = do.call("rbind", logreg_p_value_pr)
rownames(logreg_p_value_pr) = NULL
logreg_p_value_roc = do.call("rbind", logreg_p_value_roc)
rownames(logreg_p_value_roc) = NULL
logreg_p_value = rbind(logreg_p_value_pr, logreg_p_value_roc)




logreg_p_value %>% dplyr::mutate(dplyr::across(c(2,4), format, digits = 3))  %>% 
  dplyr::mutate(dplyr::across(p_value, format, digits = 4))  %>% 
  dplyr::mutate(dplyr::across(c(1,3,8), paste0, "-Weighted")) %>%  # add weighted at the end to let people know it's the logreg
  readr::write_csv(file="./sup_data/sup_data_pairwise-tis-weighted-individual.csv.gz") 



# baseline vs logreg
get_p_value_tis_basevslog = function(x, measure){
  vs = unique(x$variant_score)
  x = x[x$measure == measure, ]
  p_value = lapply(vs, function(v_score){
    auc_vs = x[x$variant_score == v_score, ]
    lr = auc_vs[auc_vs$regularization == "min",      c("auc", "id")]
    colnames(lr)[1] = c("min")
    bl = auc_vs[auc_vs$regularization == "baseline", c("auc", "id")]
    colnames(bl)[1] = c("baseline")
    a = full_join(lr, bl, by = "id")
    res = wilcox.test(a[,"min"], a[, "baseline"], paired = T)
    p_val = res$p.value
    score_median = c(median(a[,"min"]), median(a[,"baseline"]))
    names(score_median) = c("Weighted", "Mean")
    p_df = data.frame(score = v_score,
                      score1 = "Weighted", score1_median = score_median[1], 
                      score2 = "Mean", score2_median = score_median[2], 
                      disease = x$PTPE[1], 
                      measure = measure, p_value = p_val, 
                      higher_median = names(which.max(score_median)))
  })
  p_value = do.call("rbind", p_value)
  rownames(p_value) = NULL
  return(p_value)
}

comp_p_value_pr = lapply(all_SNPSNAP_logreg, get_p_value_tis_basevslog, measure = "pr")
comp_p_value_roc = lapply(all_SNPSNAP_logreg, get_p_value_tis_basevslog, measure = "roc")



get_stats2 = function(x){
  names(x) = NULL
  all = do.call("rbind", x)
  all_p = all[all$p_value < 0.05,]
  # all_p = all
  t = tapply(all_p$higher_median, all_p$score, table)
  t = do.call("rbind", t)
  t = as.data.frame(t)
  t$score = rownames(t)
  rownames(t) = NULL
  return(t)
}


c1 = get_stats2(comp_p_value_pr)
c2 = get_stats2(comp_p_value_roc)
c1$measure = "pr"
c2$measure = "roc"
cc = rbind(c1, c2)

cc$total_number = 111
cc$non_sig = 111-cc$Weighted-cc$Mean
cc2 = cc[,c(5,1,6,2,3,4)]
cc2

tb = data.frame(score = cc2$score, tissue_weighted_wins = cc2$Weighted, 
                tissue_mean_wins = cc2$Mean, tie = cc2$non_sig)
tb1 = tb[1:3,]
tb2 = tb[4:6,]
#write_csv(tb1, path="~/project/variant-scores-2/data/GWAS_catalog_project_3.0/tables/tab_tournament3.csv") 
#write_csv(tb2, path = "~/project/variant-scores-2/data/GWAS_catalog_project_3.0/tables/sup_tab_tournament3.csv")


comp_p_value_pr = do.call("rbind", comp_p_value_pr)
comp_p_value_roc = do.call("rbind", comp_p_value_roc)
comp_p_value = rbind(comp_p_value_pr, comp_p_value_roc)
rownames(comp_p_value) = NULL


comp_p_value$score1 = paste(comp_p_value$score, comp_p_value$score1, sep = "-")
comp_p_value$score2 = paste(comp_p_value$score, comp_p_value$score2, sep = "-")
comp_p_value$higher_median = paste(comp_p_value$score, comp_p_value$higher_median, sep = "-")

comp_p_value[,2:9] %>% dplyr::mutate(dplyr::across(c(2,4), format, digits = 3))  %>% 
  dplyr::mutate(dplyr::across(p_value, format, digits = 4))  %>% 
  readr::write_csv(path="./sup_data/sup_data_pairwise-tis-individual.csv.gz") 

