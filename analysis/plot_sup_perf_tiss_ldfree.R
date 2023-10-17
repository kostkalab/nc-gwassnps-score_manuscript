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

DHS_ld_free <- readRDS("./data/analysis/perf_DHS_logreg_ldfree.RDS")
DHS <- readRDS("./data/analysis/perf_DHS_logreg.RDS")

# only pick out min and baseline
DHS_ld_free = DHS_ld_free[DHS_ld_free$regularization %in% c("min", "baseline"),]
DHS         = DHS        [DHS$regularization %in% c("min", "baseline"),]


DHS$LD = "all-SNVs"
DHS_ld_free$LD = "one-SNV-per-LD"

all=rbind(DHS, DHS_ld_free)
all = all[all$measure == "pr",]

bl = all[all$regularization == "baseline" & all$run_id == 1,]
lr = all[all$regularization == "min",]

lr$PTPE = gsub("_"," ",lr$PTPE)
lr$PTPE = str_wrap(lr$PTPE,15)

bl$PTPE = gsub("_"," ",bl$PTPE)
bl$PTPE = str_wrap(bl$PTPE,15)

ph1 = unique(lr$PTPE)[1:56]
ph2 = unique(lr$PTPE)[57:111]
p1 = ggplot(lr[lr$PTPE %in% ph1,], aes(LD,auc,col=LD)) + 
  geom_quasirandom(size = 1) + 
  facet_wrap(~PTPE, ncol = 8) + 
  geom_point(data = bl[bl$PTPE %in% ph1,], aes(LD, auc), shape = 8, color = "black") +
  scale_color_brewer(palette="Dark2")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(b=12,r=5,t=7,l=5, "pt")) + 
  ggtitle(paste("all-SNVs vs one-SNV-per-LD-block")) + 
  theme(legend.position = "none")

p2 = ggplot(lr[lr$PTPE %in% ph2,], aes(LD,auc,col=LD)) + 
  geom_quasirandom(size = 1) + 
  facet_wrap("PTPE", ncol = 8) + 
  geom_point(data = bl[bl$PTPE %in% ph2,], aes(LD, auc), shape = 8, color = "black") +
  scale_color_brewer(palette="Dark2")+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(b=12,r=5,t=7,l=5, "pt")) + 
  ggtitle(paste("all-SNVs vs one-SNV-per-LD-block, continued")) + 
  theme(legend.position = "none")


pdf(file="./plot/tiss_perf_ldfree.pdf", width = 11, height = 14) # 9 and 11 is also good
p1
p2
dev.off()
