library(pROC)
library(ggplot2)
library(reshape2)
library(glmnet)
library(pamr)
library(ggbeeswarm)
library(PRROC)


source("./model/fun_logreg_chr_holdout.R")


## read in file
DHS_SNPSNAP  = readRDS("./data/anno_DHS_catalog_cutoff100_snpsnap.RDS")
chrom_heldout = readRDS("./data/catalog_cutoff100_snpsnap_chrom_heldout.RDS")

# annotate file using chrom_heldout
DHS_SNPSNAP$id = paste0(DHS_SNPSNAP$snpID, DHS_SNPSNAP$phenotype)
chrom_heldout$id = paste0(chrom_heldout$snpID, chrom_heldout$phenotype)
m = match(DHS_SNPSNAP$id, chrom_heldout$id)
DHS_SNPSNAP$chrom_heldout = chrom_heldout$chrom_heldout[m]


multiPheno = function(file, times, type){
  stopifnot(type %in% c("chrom", "random"))
  
  dat = file
  
  col         = colnames(dat)[grep("E[0-9]+", colnames(dat))]
  dat = dat[!is.na(rowSums(dat[,col])),]
  dat$phenotype = tolower(gsub(" ","_",dat$phenotype))
  
  PTPE = unique(dat$phenotype)
  if (type == "chrom"){
    multi_p = lapply(PTPE, rep_procPheno_chr, dat = dat)
  } else if (type == "random"){
    multi_p = lapply(PTPE, rep_procPheno_random, times = times, dat = dat)
  }
  
  multi_p = do.call("rbind", multi_p)
  return(multi_p)
}




## chromosome held out. Test set 1/5, train set 4/5
DHS_chr  = multiPheno(DHS_SNPSNAP, type = "chrom")
saveRDS(DHS_chr, "./data/analysis/perf_DHS_chrom.RDS")

## random held out. Test set 1/5, train set 4/5
DHS_random  = multiPheno(DHS_SNPSNAP, times = 10, type = "random")
saveRDS(DHS_random, "./data/analysis/perf_DHS_random.RDS")


DHS_cv = readRDS("~/project/variant-scores-2/analysis/GWAS_catalog_3.0/aggre/result_data/DHS_SNPSNAP_logreg.RDS")
DHS_random = readRDS("~/project/variant-scores-2/analysis/GWAS_catalog_3.0/aggre/result_data/DHS_random.RDS")
DHS_chr = readRDS("~/project/variant-scores-2/analysis/GWAS_catalog_3.0/aggre/result_data/DHS_chrom.RDS")


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
pdf(file="~/project/variant-scores-2/analysis/GWAS_catalog_3.0/aggre/plot/chromosome_held_out2.pdf", width = 15, height = 19)
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

# generate the table
DHS_df_median = t(tapply(DHS_all_pr$auc, list(DHS_all_pr$regularization, DHS_all_pr$PTPE), median))
DHS_df_sd = t(tapply(DHS_all_pr$auc, list(DHS_all_pr$regularization, DHS_all_pr$PTPE), sd))
colnames(DHS_df_sd)[4:6] = c("CV_LR(sd)", "Random_B(sd)", "Random_LR(sd)")
DHS_df = cbind(DHS_df_median, DHS_df_sd[,4:6])
DHS_df = as.data.frame(DHS_df)
DHS_df$disease = rownames(DHS_df)
DHS_df = DHS_df[,c(10,1,2,3,4,7,5,8,6,9)]

readr::write_csv(DHS_df, file = "~/project/variant-scores-2/data/GWAS_catalog_project_3.0/tables/sup_data_perf-chrom-heldout.csv.gz")
