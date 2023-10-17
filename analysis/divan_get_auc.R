library(PRROC)
library(pROC)

source("./model/fun_logreg_divan.R")

catalog_DIVAN = readRDS("./data/anno_divan_DIVAN.RDS")
catalog_DHS   = readRDS("./data/anno_divan_DHS.RDS")
catalog_org   = readRDS("./data/anno_divan_GenoCanyon.RDS")
test_boostrap = readRDS("./data/analysis/bootstrap_samples_divan29.RDS")

pheno = tolower(gsub(" ","_",names(test_boostrap)))
names(test_boostrap) = pheno

# give catalog DIVAN names
for (i in seq_along(catalog_DIVAN)){
  names(catalog_DIVAN)[i] = tolower(gsub(" ","_",catalog_DIVAN[[i]]$trait))[1]
}
catalog_DIVAN = catalog_DIVAN[names(catalog_DIVAN) %in% pheno]
catalog_DHS$phenotype = tolower(gsub(" ","_",catalog_DHS$phenotype))
catalog_DHS = catalog_DHS[catalog_DHS$phenotype %in% pheno,]
catalog_org$trait = tolower(gsub(" ","_",catalog_org$trait))
catalog_org = catalog_org[catalog_org$trait %in% pheno]



DHS_auc = lapply(pheno, procPhenotype_oldnew, dat = catalog_DHS)
auc_DHS = do.call("rbind", DHS_auc)
auc_DHS$measure = "min"
rownames(auc_DHS) = NULL


org_auc = function(x){
  trait = tolower(gsub(" ","_",as.character(x$trait[1])))
  print(trait)
  X = x[x$divan_new == 1]
  X = X[!is.na(X$score.mean)] # remove NA
  y = X$type
  bts = data.frame(test_boostrap[[trait]], stringsAsFactors = F)
  res = sapply(bts, function(x){
    X_samp = X[X$snpID %in% x,]
    auroc = roc(X_samp$type, X_samp$score.mean)$auc
    aupr = pr.curve(scores.class0=X_samp$score.mean,weights.class0=X_samp$type)$auc.integral
    res1 = c(auroc, aupr)
    names(res1) = c("roc", "pr")
    return(res1)
  })
  aucd = melt(res)
  colnames(aucd) = c("type", "id", "auc")
  aucd$PTPE = trait
  return(aucd)
}


auc_divan = lapply(catalog_DIVAN, org_auc)
auc_divan = do.call("rbind", auc_divan)
auc_divan$measure = "DIVAN"
rownames(auc_divan) = NULL

# org performance
colnames(mcols(catalog_org))[5] = "score.mean" # run this line with careful, sometimes it's not [8]!!!!!!
catalog_org = split(catalog_org, catalog_org$trait)

auc_org = lapply(catalog_org, org_auc)
auc_org = do.call("rbind", auc_org)
auc_org$measure = "GenoCanyon"
rownames(auc_org) = NULL

auc_all = rbind(auc_divan, auc_org, auc_DHS)

saveRDS(auc_all, "./data/analysis/perf_all_divan.RDS")

