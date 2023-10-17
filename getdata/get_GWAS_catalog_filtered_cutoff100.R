library(rols)
library(tibble)
library(stringr)
library(dplyr)
library(tidyr)

library(ontologyIndex) # remember to use version 2.2 to have the get_OWL function
file = "./external_data/efo_ontology/3.24.0/efo_3_24_0.obo"
efo = get_ontology(file, propagate_relationships = "is_a", extract_tags = "minimal")


catalog_aggre = readRDS("./data/GWAS_catalog_noncoding_hg19_uniquetrait_propagate.RDS")
catalog       = readRDS("./data/GWAS_catalog_noncoding_hg19_uniquetrait.RDS")
catalog_efo   = catalog[grep("EFO", catalog$UNIQUE.URI)] # 194759 -> 184785 SNPs  # 236675 -> 225698 # 236470 -> 225506 # keep only efo terms

# keep only unique SNPs/phenotype
catalog = tapply(catalog_efo, catalog_efo$UNIQUE.URI, function(x){
  x1 = x[!duplicated(x)] # don't keep duplicated SNPs per phenotype
  return(x1)
}) 


catalog_100_aggre = catalog_aggre[sapply(catalog_aggre, length) >= 100]
catalog_100       = catalog[sapply(catalog, length) >= 100]

nms = names(catalog_100_aggre)


d_trm_efo = sapply(nms, function(x){
  # print(x)
  anc = get_ancestors(efo, x)
  print(is.null(anc))
  if(!is.null(anc)){
    re = "EFO:0000408" %in% anc
  }else{
    re = T
  }
  
  return(re)
})

nms = nms[d_trm_efo]

catalog_100_d       = catalog_100[names(catalog_100) %in% nms]
catalog_100_aggre_d = catalog_100_aggre[d_trm_efo]

saveRDS(catalog_100_aggre_d, "./data/GWAS_catalog_filtered_cutoff100.RDS")
saveRDS(catalog_100_d,       "./data/GWAS_catalog_filtered_cutoff100_before_propagation.RDS")

# generate a file that contains unique snps to match using snpsnap
cata_unlist = catalog_100_aggre_d
names(cata_unlist) = NULL
cata_unlist = do.call("c", cata_unlist)
chr = gsub("chr","",seqnames(cata_unlist))
chr = gsub("X", "23", chr)

g.coord = paste0(chr,":",start(cata_unlist))
g.coord = c("lead_snp",unique(g.coord))

# delete y chromosome
g.coord = g.coord[!grepl("Y", g.coord)]
write.table(g.coord, 
            file = "~/snpsnap/to_dennis_v2_2020_04_27/GWAS_inputsnps_all_ct100.txt", 
            sep = "\n", row.names = F, col.names = F, quote = F)

