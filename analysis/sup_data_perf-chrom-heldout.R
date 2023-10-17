library(reshape2)


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



# generate the table
DHS_df_median = t(tapply(DHS_all_pr$auc, list(DHS_all_pr$regularization, DHS_all_pr$PTPE), median))
DHS_df_sd = t(tapply(DHS_all_pr$auc, list(DHS_all_pr$regularization, DHS_all_pr$PTPE), sd))
colnames(DHS_df_sd)[4:6] = c("CV_LR(sd)", "Random_B(sd)", "Random_LR(sd)")
DHS_df = cbind(DHS_df_median, DHS_df_sd[,4:6])
DHS_df = as.data.frame(DHS_df)
DHS_df$disease = rownames(DHS_df)
DHS_df = DHS_df[,c(10,1,2,3,4,7,5,8,6,9)]

readr::write_csv(DHS_df, file = "./sup_data/sup_data_perf-chrom-heldout.csv.gz")
