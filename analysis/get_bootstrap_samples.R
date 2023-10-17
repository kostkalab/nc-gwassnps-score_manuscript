library(GenomicRanges)
library(sampler)

catalog_org  = readRDS("./data/catalog_cutoff100_snpsnap.RDS")


disease = catalog_org[catalog_org$type == 1]
sc = c("CADD", "eigen", "GenoCanyon", "GWAVA", "LINSIGHT")
t = table(catalog_org[catalog_org$type == 1]$phenotype)
phenotypes = names(t[t>=100])


bootstrap_samples = lapply(phenotypes, function(x){
  set.seed(0525)
  print(x)
  catalog = data.frame(mcols(catalog_org[catalog_org$phenotype == x]))
  get_sub = function(x){
    sub_cat = sampler::ssamp(df=x, n=floor(nrow(x) * 0.9), strata=type) # sample
    snpID   = sub_cat$snpID
    return(snpID)
  }
  sub_30 = replicate(30, get_sub(catalog))
  colnames(sub_30) = 1:30
  return(sub_30)
})

names(bootstrap_samples) = phenotypes

saveRDS(bootstrap_samples, "./data/analysis/bootstrap_samples_30_0.9.RDS")
