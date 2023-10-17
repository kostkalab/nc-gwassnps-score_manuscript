library(rols)
library(tibble)
library(stringr)
library(dplyr)
library(tidyr)

library(ontologyIndex) # remember to use version 2.2 to have the get_OWL function
file = "~/project/variant-scores_manuscript/external_data/efo_ontology/3.24.0/efo_3_24_0.obo"
efo = get_ontology(file, propagate_relationships = "is_a", extract_tags = "minimal") # i'm pretty sure rols package use is_a not part_of

catalog_all = readRDS("~/project/variant-scores_manuscript/data/GWAS_catalog_noncoding_hg19_uniquetrait.RDS")

# keep only efo trait
catalog_efo = catalog_all[grep("EFO", catalog_all$UNIQUE.URI)] # 236470 -> 225506


# keep only unique SNPs/phenotype
catalog_unique = tapply(catalog_efo, catalog_efo$UNIQUE.URI, function(x){
  x1 = x[!duplicated(x)] # don't keep duplicated SNPs per phenotype
  return(x1)
}) # can use this directly


# go over each trait and add aggregate
# the following function uses efo.obo and ontologyIndex package to get the descendants. 
# this is also faster than rols package. 
catalog_aggre = lapply(catalog_unique, function(x){
  AGGRE.URI   = x$UNIQUE.URI[1]
  print(AGGRE.URI)
  AGGRE.TRAIT = x$UNIQUE.TRAIT[1]
  # trm      = term("efo", AGGRE.URI)
  # desc     = descendants2(trm)
  desc = get_descendants(efo, AGGRE.URI, exclude_roots = T)
  print(length(desc))
  x$aggre = F
  if (length(desc) > 0){
    c_desc   = catalog_efo[catalog_efo$UNIQUE.URI %in% desc]
    if (length(c_desc) > 0){
      c_desc$aggre = T
      c_aggre  = c(x, c_desc)
      c_aggre  = c_aggre[!duplicated(c_aggre)] # delete duplicated SNPs after aggregation (delete aggregated ones)
    } else {c_aggre = x }
    
  } else { c_aggre = x}
  
  c_aggre$AGGRE.URI = AGGRE.URI
  c_aggre$AGGRE.TRAIT = AGGRE.TRAIT
  
  return(c_aggre)
})


saveRDS(catalog_aggre, "~/project/variant-scores_manuscript/data/GWAS_catalog_noncoding_hg19_uniquetrait_propagate.RDS")

#### compared with the old ones i used rols, it is exactly the same!! hooray!! 
