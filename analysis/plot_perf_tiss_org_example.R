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


# pick the organism level score to use
org_score = "GenoCanyon"

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


## read in organism level scores
org_auc <- readRDS("./data/analysis/perf_org.RDS")
org_auc$PTPE = firstup(tolower(org_auc$PTPE))
org_snpsnap = org_auc[org_auc$match == "SNPsnap",]
org_snpsnap = org_snpsnap[,c("variant_score", "measure", "match", "auc", "PTPE")]
colnames(org_snpsnap) = c("regularization", "measure", "run_id", "auc", "PTPE")


plot_tis_org = function(df, org, name, measure){
  ## baseline_tis performance
  bl = df[df$regularization == "baseline" & df$run_id == 1 & df$measure == measure & df$variant_score == "DHS",]
  bl = bl[order(bl$PTPE),]
  bl = bl[order(bl$variant_score),]
  
  ## org performance
  org_1 = org[org$regularization == org_score & 
                org$measure == measure,]
  
  ## logreg_median performance
  ml = df[df$regularization == name & df$measure == measure & df$variant_score == "DHS",]
  ml_median = tapply(ml$auc, list(ml$PTPE, ml$variant_score), median)
  ml_median_long = reshape2::melt(ml_median)
  colnames(ml_median_long) = c("PTPE", "variant_score", "auc")
  
  ## get the vector to sort ml performance
  DHS_median = ml_median[,"DHS"]
  diff = (DHS_median - org_1$auc)[order(DHS_median - org_1$auc, decreasing = T)]
  
  v1 = diff[diff >= 0.05][1:4] # total 24
  v2 = diff[diff >= 0 & diff < 0.05][c(1:4,40:43)] # total 74
  v3 = diff[diff < 0][10:13] # total 13
  
  df_name = enframe(c(v1,v2,v3))
  df_name$group = rep(1:4, each = 4)
  df_name$group = rep(c("Best Improvement", "Medium Improvement", "Comparable Performance", paste(org_score, "Better")), each = 4)
  df_name$group = factor(df_name$group, levels = c("Best Improvement", "Medium Improvement", "Comparable Performance", paste(org_score, "Better")))
  
  DHS_sort = df_name$name
  ml2 = ml[ml$PTPE %in% df_name$name,]
  ml2$PTPE = factor(ml2$PTPE, levels = DHS_sort)
  
  # df takes in tissue specific scores
  # or takes in organism score
  or = org[org$measure == measure, ]
  scores = unique(or$regularization)
  orr = or[or$regularization == org_score, c("auc", "PTPE")]
  orr3 = orr[orr$PTPE %in% DHS_sort,]
  orr3$variant_score = "DHS"
  
  colnames(df_name)[1] = "PTPE"
  orr4 = full_join(orr3, df_name[,c("PTPE", "group")])
  ml3 = full_join(ml2, df_name[,c("PTPE", "group")])
  ml3$PTPE = factor(ml3$PTPE, levels = DHS_sort)
  
  # combine ml3 and orr3 together
  ml3 = ml3[,c("auc", "PTPE","variant_score","group")]
  orr4$variant_score = "GenoCanyon"
  
  ml4 = rbind(ml3, orr4)
  cols <- c("DHS"="#f04546","GenoCanyon"="#3591d1")
  p   = ggplot(ml4, aes(PTPE,auc,group=variant_score)) + 
    geom_quasirandom(aes(shape=variant_score, color=variant_score), size = 1) + 
    scale_shape_manual(name = "Variant score",values=c(16, 5))+
    scale_color_manual(name = "Variant score",values=c('#F8766D','#000000'))+ 
    # scale_size_manual(values=c('#F8766D','#000000'))+ 
    # geom_point(data = orr4, aes(PTPE, auc), col = "black", shape = 5, size = 0.7) + 
    facet_wrap(c("group"), nrow = 1, scales = "free_x") + 
    theme_bw()+ 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          plot.margin = margin(0.2, 0.2, 0.2, 0.5, "cm")) + 
    ggtitle("") + 
    theme(axis.text=element_text(size=8),
          strip.background = element_blank()) + 
    ylab("Average precision")+
    xlab("Disease term") 
  return(p)
}

p2 = plot_tis_org(df = all_p_SNPSNAP, org = org_snpsnap, name = "min", measure = "pr")

pdf("./plot/perf-tissue-vs-organism-individual-sub.pdf", width = 8, height = 3.7)
p2
dev.off()
