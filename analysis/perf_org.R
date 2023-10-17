library(PRROC)
library(pROC)
library(GenomicRanges)



source("./analysis/fun_perf_org.R")

tss <- readRDS("./data/anno_org_catalog_cutoff100_tss.RDS")
snpsnap_tss <- readRDS("./data/anno_org_catalog_cutoff100_snpsnaptss.RDS")
snpsnap <- readRDS("./data/anno_org_catalog_cutoff100_snpsnap.RDS")
random <- readRDS("./data/anno_org_catalog_cutoff100_random.RDS")


get_auc = function(catalog_org, matched_method){
  case    = catalog_org[catalog_org$type == 1]
  control = catalog_org[catalog_org$type == 0] 
  
  sc = c("CADD", "eigen", "GenoCanyon", "GWAVA", "LINSIGHT")
  t = table(catalog_org$phenotype)
  phenotypes = names(t[t>30]) # why 30? 
  
  org_auc = lapply(phenotypes, make_roc_org, ctrl = control, case = case, sc = sc) # why cannot use GenoCanyon? 
  org_auc = do.call("rbind", org_auc)
  org_auc$match = matched_method
  return(org_auc)
}

a1 = get_auc(tss, "TSS")
a2 = get_auc(snpsnap_tss, "SNPsnap_TSS")
a3 = get_auc(snpsnap, "SNPsnap")
a4 = get_auc(random, "random")

org_auc = rbind(a1, a2, a3, a4)

saveRDS(org_auc, "./data/analysis/perf_org.RDS")
