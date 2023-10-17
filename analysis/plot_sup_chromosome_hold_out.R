
library(ggplot2)
library(reshape2)

library(pamr)
library(ggbeeswarm)



DHS_cv = readRDS("./data/analysis/perf_DHS_logreg.RDS")
DHS_random = readRDS("./data/analysis/perf_DHS_random.RDS")
DHS_chr = readRDS("./data/analysis/perf_DHS_chrom.RDS")

# organize so that i can combine them
DHS_random = DHS_random[DHS_random$regularization %in% c("min", "baseline"),]
DHS_random$regularization = as.character(DHS_random$regularization)
DHS_random[DHS_random$regularization == "baseline", ]$regularization = "Random_B"
DHS_random[DHS_random$regularization == "min", ]$regularization = "Random_LR"

DHS_cv = DHS_cv[DHS_cv$regularization %in% c("min", "baseline"),]
DHS_cv$regularization = as.character(DHS_cv$regularization)
DHS_cv[DHS_cv$regularization == "baseline", ]$regularization = "CV_B"
DHS_cv[DHS_cv$regularization == "min", ]$regularization = "CV_LR"
DHS_cv = DHS_cv[,1:5]

DHS_chr = DHS_chr[DHS_chr$regularization %in% c("min", "baseline"),]
DHS_chr$regularization = as.character(DHS_chr$regularization)
DHS_chr[DHS_chr$regularization == "baseline", ]$regularization = "Chr_B"
DHS_chr[DHS_chr$regularization == "min", ]$regularization = "Chr_LR"
DHS_chr$run_id = 1

# combine them 
DHS_all = rbind(DHS_cv, DHS_random, DHS_chr)
DHS_all_pr = DHS_all[DHS_all$measure == "pr",] # get the PR 

# make them into two groups to plot in two pages
diseases = unique(DHS_chr$PTPE)
PTPE1 = diseases[1:56]
PTPE2 = diseases[57:111]
am1 = DHS_all_pr[DHS_all_pr$PTPE %in% PTPE1,]
am2 = DHS_all_pr[DHS_all_pr$PTPE %in% PTPE2,]

# color = brewer.pal(12, "Paired")[c(1:2, 5:6, 9:10)] # if i want to change color i can modify this line
pdf(file="./plot/chromosome_held_out.pdf", width = 15, height = 19)
p = ggplot(am1,       aes(regularization,auc,col = regularization)) + geom_boxplot()+ 
  geom_quasirandom() + 
  geom_hline(yintercept=1/11, linetype="dashed") + 
  facet_wrap("PTPE", ncol = 7) + 
  scale_color_brewer(palette="Paired") + 
  # scale_color_manual(values=color) +  # modify this line if i want to change color
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("Area Under Precision Recall Curve") + 
  xlab("held out method")
p
p = ggplot(am2,       aes(regularization,auc,col = regularization)) + geom_boxplot()+ 
  geom_quasirandom() + 
  geom_hline(yintercept=1/11, linetype="dashed") + 
  facet_wrap("PTPE", ncol = 7) + 
  scale_color_brewer(palette="Paired") + 
  # scale_color_manual(values=color) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("Area Under Precision Recall Curve") + 
  xlab("held out method")
p
dev.off()

