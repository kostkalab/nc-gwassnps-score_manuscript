library(pROC)
library(ggplot2)
library(reshape2)
library(glmnet)
library(pamr)
library(ggbeeswarm)
library(PRROC)
library(dplyr)
library(tidyr)
library(kableExtra)
library(gridExtra)
library(tibble)

## read in data and manage it
# read in auc for logistic regression and xgbtree
DHS_SNPSNAP_logreg <- readRDS("./data/analysis/perf_DHS_logreg.RDS")
Fit_SNPSNAP_logreg <- readRDS("./data/analysis/perf_Fit_logreg.RDS")
geno_SNPSNAP_logreg <- readRDS("./data/analysis/perf_geno_logreg.RDS")

all_p_SNPSNAP = rbind(DHS_SNPSNAP_logreg, Fit_SNPSNAP_logreg, geno_SNPSNAP_logreg)


## read in organism level scores
org_auc <- readRDS("./data/analysis/perf_org.RDS")

# because there are three matching strategies now, have to pick out snpsnap one. note on: 12/22/2020
org_auc$PTPE = tolower(gsub(" ","_",org_auc$PTPE))
org_snpsnap = org_auc[org_auc$match == "SNPsnap",]
org_snpsnap = org_snpsnap[,c("variant_score", "measure", "match", "auc", "PTPE")]
colnames(org_snpsnap) = c("regularization", "measure", "run_id", "auc", "PTPE")


## get the factor that sort by organism score and tissue scores difference
get_quantile = function(p, df, mdname){
  ml = df[df$regularization == mdname,]
  md = tapply(ml$auc, INDEX = list(ml$measure, ml$PTPE, ml$variant_score), FUN = quantile, probs = p)
  md_long = melt(md,varnames = c("measure", "PTPE", "variant_score"), value.name = "auc")
  md_long$regularization = p
  return(md_long)
}



qt = c(0.25, 0.5, 0.75)
logreg_SNPSNAP_quantile = lapply(qt, get_quantile, df = all_p_SNPSNAP, mdname = "min")
logreg_SNPSNAP_quantile = do.call("rbind", logreg_SNPSNAP_quantile)
logreg_SNPSNAP_quantile$regularization = paste0(logreg_SNPSNAP_quantile$regularization, "_log")



plot_tis_org = function(df, org, name, measure){
  ## baseline_tis performance
  bl = df[df$regularization == "baseline" & df$run_id == 1 & df$measure == measure,]
  bl = bl[order(bl$PTPE),]
  bl = bl[order(bl$variant_score),]
  
  ## org performance
  org_1 = org[org$regularization == "GenoCanyon" & 
                org$measure == measure,]
  
  ## logreg_median performance
  ml = df[df$regularization == name & df$measure == measure,]
  ml_median = tapply(ml$auc, list(ml$PTPE, ml$variant_score), median)
  ml_median_long = melt(ml_median)
  colnames(ml_median_long) = c("PTPE", "variant_score", "auc")
  
  ## get the vector to sort ml performance
  DHS_median = ml_median[,"DHS"]
  DHS_sort = names((DHS_median - org_1$auc)[order(DHS_median - org_1$auc, decreasing = T)])
  ml$PTPE = factor(ml$PTPE, levels = DHS_sort)
  
  # get a dataframe to plot
  ml_median_long$regularization = "logreg"
  df_all = rbind(ml_median_long, bl[,c("PTPE", "variant_score", "auc", "regularization")])
  names(org_1)[4] = "org_auc"
  b = full_join(df_all, org_1[,c("org_auc", "PTPE")], by = "PTPE")
  b$diff = b$auc - b$org_auc
  
  
  g1 <- ggplot(b, aes(y=diff,x = regularization, col = variant_score)) + 
    geom_quasirandom(size = 1) + 
    facet_wrap(c("variant_score", "regularization"), nrow = 3, scales = "free_x") +
    ylab("tissue scores - GenoCanyon") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    ggtitle("   ") +
    theme(legend.position = "none")+ 
    theme(strip.text.x = element_text(size = 6))
  
  
  # have a color vector that indicates whether it is overlap with 0.25-0.75
  org_1
  tiss_25 = logreg_SNPSNAP_quantile[logreg_SNPSNAP_quantile$regularization == "0.25_log" &
                                      logreg_SNPSNAP_quantile$measure == measure,]
  tiss_75 = logreg_SNPSNAP_quantile[logreg_SNPSNAP_quantile$regularization == "0.75_log" &
                                      logreg_SNPSNAP_quantile$measure == measure,]
  compare_df = data.frame(PTPE = tiss_25$PTPE, q_25 = tiss_25$auc, q_75 = tiss_75$auc, org = org_1$org_auc, 
                          variant_score = tiss_25$variant_score)
  compare_df$col = "black"
  compare_df$col[compare_df$q_25 >= compare_df$org] = "red"
  compare_df$col[compare_df$q_75 <= compare_df$org] = "blue"
  # i don't need to sort this by rank before, they are going to automatically do it
  
  
  # df takes in tissue specific scores
  # or takes in organism score
  or = org[org$measure == measure, ]
  scores = unique(or$regularization)
  orr = or[or$regularization == "GenoCanyon", c("auc", "PTPE")]
  orr3 = rbind(orr, orr, orr)
  n = length(unique(df$PTPE))
  orr3$variant_score = c(rep("DHS", n), rep("Fitcons2", n), rep("Genoskyline", n))
  p   = ggplot(ml, aes(PTPE,auc,col=variant_score)) + 
    geom_quasirandom(size = 1) + 
    geom_point(data = orr3, aes(PTPE, auc), col = compare_df$col, shape = 5, size = 0.7) + 
    facet_wrap("variant_score", nrow = 3) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          plot.margin = margin(0.2, 0.2, 0.2, 1.5, "cm")) + 
    ggtitle(paste0("Tissue-weighted-DHS vs. GenoCanyon ", measure)) + 
    theme(legend.position = "none", axis.text=element_text(size=8)) + 
    ylab(paste0("AU", toupper(measure))) 
  both = list(g1, p)
  return(both)
}

p1 = plot_tis_org(df = all_p_SNPSNAP, org = org_snpsnap, name = "min", measure = "roc")
p2 = plot_tis_org(df = all_p_SNPSNAP, org = org_snpsnap, name = "min", measure = "pr")


pdf("./plot/tiss-vs-org-perf-individual-pr.pdf", width = 13, height = 7)
# p1[[2]]
p2[[2]]
dev.off()