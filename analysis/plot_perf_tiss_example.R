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

# read in SNPSNAP tissue specific scores
DHS_SNPSNAP_logreg <- readRDS("./data/analysis/perf_DHS_logreg.RDS")
Fit_SNPSNAP_logreg <- readRDS("./data/analysis/perf_Fit_logreg.RDS")
geno_SNPSNAP_logreg <- readRDS("./data/analysis/perf_geno_logreg.RDS")


all_p_SNPSNAP = rbind(DHS_SNPSNAP_logreg, Fit_SNPSNAP_logreg, geno_SNPSNAP_logreg)

# change disease names into the format that we want: 
## remove underscore but spaces
## capitalize the first letter

firstup = function(x){
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
all_p_SNPSNAP$PTPE = firstup(gsub("_", " ",all_p_SNPSNAP$PTPE))

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
  diff = (DHS_median - DHS_baseline)[order(DHS_median - DHS_baseline, decreasing = T)]
  
  v1 = diff[diff >= 0.10][1:4] # length 8
  
  # v2 = diff[diff >= 0.05 & diff < 0.10][7:10] # length 27
  v2 = v2 = diff[diff >= 0.05 & diff < 0.10][12:15]
  # v3 = diff[diff < 0.05][30:33]
  v3 = diff[diff < 0.05][73:76] # length 76
  
  df_name = enframe(c(v1,v2,v3))
  df_name$group = rep(c("Best Improvement", "Middle Improvements", "Least Improvements"), each = 4)
  df_name$group = factor(df_name$group, levels = c("Best Improvement", "Middle Improvements", 
                                                   "Least Improvements"))
  
  DHS_sort = df_name$name
  ml2 = ml[ml$PTPE %in% df_name$name,]
  ml2$PTPE = factor(ml2$PTPE, levels = DHS_sort)
  bl2 = bl[bl$PTPE %in% DHS_sort,]
  colnames(df_name)[1] = "PTPE"
  ml3 = full_join(ml2, df_name[,c("PTPE"  ,"group")], by = "PTPE")
  bl3 = full_join(bl2, df_name[,c("PTPE"  ,"group")], by = "PTPE")
  ml3$PTPE = factor(ml3$PTPE, levels = DHS_sort)
  p = ggplot(ml3, aes(PTPE,auc,col=variant_score)) + 
    geom_quasirandom(size = 1) + 
    # facet_wrap(.~ variant_score + group, nrow = 3, scales = "free_x") + 
    facet_grid(rows = vars(variant_score),cols = vars(group),scales = "free_x",space = "free") +
    geom_point(data = bl3, aes(PTPE, auc), col = "black", shape = 8) +
    ggtitle("") + 
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(b=1,r=1,t=1,l=1, "pt"),
          strip.background = element_blank()) + 
    theme(legend.position = "none") + 
    ylab("Average precision") + 
    xlab("Disease term")
  
  return(p)
}



p2 = plot_45(df = all_p_SNPSNAP, name = "min", measure = "pr")
plot_45(DHS_SNPSNAP_logreg, "min", "pr")

pdf(file="./plot/perf-tissue-individual-sub.pdf", width = 6, height = 6)
p2
dev.off()


# only plot DHS - use this for the ASHG poster, not in manuscript
pdf(file="./plot/ASHG_only_DHS.pdf", width = 6, height = 3.6)
plot_45(DHS_SNPSNAP_logreg, "min", "pr")
dev.off()

