library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(data.table)

org_auc <- readRDS("./data/analysis/perf_org.RDS")
org_auc = org_auc[org_auc$match == "SNPsnap",]
org_auc$PTPE = as.character(org_auc$PTPE)

line = data.frame(measure = c("roc", "pr"), Z = c(0.5, 1/11))
org_auc <- data.table(org_auc)
org_auc[measure == "pr",y_min := 0.05]
org_auc[measure == "pr",y_max := 0.30]
org_auc[measure == "roc",y_min := 0.47]
org_auc[measure == "roc",y_max := 0.76]



# plot snpsnap in 45 diseases
# png(file="~/projects/variant-scores-2/analysis/results/organism_performance_45_ROC.png", width = 440, height = 350)
auc_labeller = c(
  pr = "Precision Recall Curve",
  roc = "Receiver Operating Characteristic Curve"
)

pdf(file="./plot/perf-organism-level.pdf", width = 8, height = 5)
p = ggplot(org_auc[org_auc$match == "SNPsnap",],       aes(variant_score,auc,col = variant_score)) + geom_boxplot() +
  # ggtitle("Organism-level score performance") + 
  geom_quasirandom(dodge.width=0.75, size = 0.8) + 
  facet_wrap("measure", scales = "free_y",
             labeller = as_labeller(auc_labeller)) + 
  geom_hline(data = line, aes(yintercept = Z), linetype="dashed") + 
  geom_blank(aes(y = y_min)) + 
  geom_blank(aes(y = y_max)) + 
  theme_bw() + 
  xlab("") + 
  ylab("AUC") + 
  theme(legend.position = "none", 
        strip.background = element_blank())
p
dev.off()


## # only plot precision recall curve
## pdf(file="~/projects/variant-scores-2/analysis/GWAS_catalog_3.0/aggre/plot/org3.pdf", width = 6, height = 6)
## p = ggplot(org_auc[org_auc$match == "SNPsnap" & org_auc$measure == "pr",],       aes(variant_score,auc,col = variant_score)) + geom_boxplot() +
##   ggtitle("Organism-level score performance") + geom_quasirandom() + 
##   geom_hline(yintercept=1/11, linetype="dashed") + 
##   geom_blank(aes(y = y_min)) + 
##   geom_blank(aes(y = y_max)) + 
##   ylab("Area Under Precision Recall Curve")
## p
## dev.off()
## 
## # proposal: write the average auc
## t = tapply(org_auc$auc, list(org_auc$variant_score, org_auc$measure), mean)
## t2 = tapply(org_auc$auc, list(org_auc$variant_score, org_auc$measure), median)## 