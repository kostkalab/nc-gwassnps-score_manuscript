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


## read in data and manage it
# read in auc for logistic regression and xgbtree

DHS_SNPSNAP_logreg <- readRDS("./data/analysis/perf_DHS_logreg.RDS")
Fit_SNPSNAP_logreg <- readRDS("./data/analysis/perf_Fit_logreg.RDS")
geno_SNPSNAP_logreg <- readRDS("./data/analysis/perf_geno_logreg.RDS")


all_p_SNPSNAP = rbind(DHS_SNPSNAP_logreg, Fit_SNPSNAP_logreg, geno_SNPSNAP_logreg)


# sort the same way like phegen

plot_45 = function(df, name, measure){
  
  ## regularly start from here
  bl = df[df$regularization == "baseline" & df$run_id == 1 & df$measure == measure,]
  ml = df[df$regularization == name & df$measure == measure,]
  bl = bl[order(bl$PTPE),]
  bl = bl[order(bl$variant_score),]
  
  # sort the plot by DHS_logreg_median - DHS_baseline
  DHS_baseline = bl[bl$variant_score == "DHS",]
  DHS_baseline = DHS_baseline[order(DHS_baseline$PTPE),"auc"]
  ml_median = tapply(ml$auc, list(ml$PTPE, ml$variant_score), median)
  DHS_median = ml_median[,"DHS"]
  DHS_sort = names((DHS_median - DHS_baseline)[order(DHS_median - DHS_baseline, decreasing = T)])
  ml$PTPE = factor(ml$PTPE, levels = DHS_sort)
  
  p = ggplot(ml, aes(PTPE,auc,col=variant_score)) + 
    geom_quasirandom(size = 1) + 
    facet_wrap("variant_score", nrow = 3) + 
    geom_point(data = bl, aes(PTPE, auc), col = "black", shape = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(b=12,r=1,t=7,l=50, "pt")) + 
    ggtitle(paste("Tissue-weighted vs. Tissue-mean", measure)) + 
    theme(legend.position = "none")
  
  ml_median_long = melt(ml_median)
  colnames(ml_median_long) = c("PTPE", "variant_score", "auc")
  df_all = data.frame(diff = ml_median_long$auc - bl$auc, variant_score = bl$variant_score, PTPE = bl$PTPE, 
                      type = "baseline")
  g1 <- ggplot(df_all, aes(y=diff,x = type, col = variant_score)) + 
    geom_quasirandom(size = 1) + 
    facet_wrap(c("variant_score"), nrow = 3, scales = "free_x") +
    ylab("tissue score logreg median - baseline") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    ggtitle("   ") +
    theme(legend.position = "none")+ 
    theme(strip.text.x = element_text(size = 6))
  
  both = list(g1, p)
  
  return(both)
}


# p1 = plot_45(df = all_p_SNPSNAP, name = "min", measure = "roc")
p2 = plot_45(df = all_p_SNPSNAP, name = "min", measure = "pr")


pdf(file="./plot/tiss-perf-individual-pr.pdf", width = 15, height = 7)
# p1[[2]]
p2[[2]]
dev.off()


