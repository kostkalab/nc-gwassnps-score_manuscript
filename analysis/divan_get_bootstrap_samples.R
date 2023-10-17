library(GenomicRanges)
library(ontologyIndex)


catalog_DHS   = readRDS("./data/anno_divan_DHS.RDS")
mesh_efo_mapping_45 <- readRDS("./data/mesh_efo_mapping_45.RDS")


file = "./external_data/efo_ontology/3.24.0/efo_3_24_0.obo"
efo = get_ontology(file, propagate_relationships = "is_a", extract_tags = "minimal") # i'm pretty sure rols package use is_a not part_of

disease_desc = get_descendants(efo, "EFO:0000408")
non_disease_trait = mesh_efo_mapping_45[!mesh_efo_mapping_45$mapped_label %in% efo$name[disease_desc],]


# some catalog_DHS have very few snps
# only keep phenotypes that have >= 10 new snps, and >= 50 old snps
tmp = catalog_DHS[catalog_DHS$type == 1,]
t = tapply(as.factor(tmp$divan_new), tmp$phenotype, table)
t = do.call("rbind", t)
t = t[t[,"0"] >= 50 & t[,"1"] > 20, ]
colnames(t) = c("divan_old", "divan_new")
pheno = rownames(t)

# delete non-disease traits.
pheno = pheno[!pheno %in% non_disease_trait$label]
catalog_test = data.frame(catalog_DHS[catalog_DHS$divan_new == 1, c("phenotype", "divan_new", "type","snpID")])

bootstrap_samples = lapply(pheno, function(x){
  set.seed(0525)
  print(x)
  catalog_p = catalog_test[catalog_test$phenotype == x,]
  get_sub = function(x){
    sub_cat = sampler::ssamp(df=x, n=floor(nrow(x) * 0.9), strata=type) # sample
    snpID   = sub_cat$snpID
    return(snpID)
  }
  sub_30 = replicate(30, get_sub(catalog_p))
  colnames(sub_30) = 1:30
  return(sub_30)
})

names(bootstrap_samples) = pheno
saveRDS(bootstrap_samples, "./data/analysis/bootstrap_samples_divan29.RDS")
