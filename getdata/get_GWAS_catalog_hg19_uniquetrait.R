library(rols)
library(tibble)
library(stringr)
library(dplyr)
library(tidyr)

library(ontologyIndex) # remember to use version 2.2 to have the get_OWL function
file = "./external_data/efo_ontology/3.24.0/efo_3_24_0.obo"
efo = get_ontology(file, propagate_relationships = "is_a", extract_tags = "minimal") # i'm pretty sure rols package use is_a not part_of

GWAS_catalog_noncoding_hg19 <- readRDS("./data/GWAS_catalog_noncoding_hg19.RDS")
GWAS_catalog_noncoding_hg19$MAPPED_TRAIT

## make this into a tibble object
GWAS_catalog = data.frame(mcols(GWAS_catalog_noncoding_hg19), seqnames = seqnames(GWAS_catalog_noncoding_hg19), 
                          start = start(GWAS_catalog_noncoding_hg19), 
                          end = end(GWAS_catalog_noncoding_hg19))
GWAS_catalog = as_tibble(GWAS_catalog)


# find out the trait name that contains , in itself
s1 = str_count(GWAS_catalog$MAPPED_TRAIT, ", ")
s2 = str_count(GWAS_catalog$MAPPED_TRAIT_URI, ", ")

# split the all dataset into two part
# part 1: the trait that has , in itself
catalog_1 = GWAS_catalog[s1 - s2 > 0, ]
catalog_1 = catalog_1 %>%
  rowwise() %>%
  mutate(MAPPED_TRAIT_URI = str_split(MAPPED_TRAIT_URI,pattern = ", ")) %>%
  unnest(c(MAPPED_TRAIT_URI)) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
catalog_1$MAPPED_TRAIT_URI = basename(catalog_1$MAPPED_TRAIT_URI)
catalog_1$MAPPED_TRAIT_URI = gsub("_", ":", catalog_1$MAPPED_TRAIT_URI)

# use rols package to name it to replace the original one. 
s = split(catalog_1, f = catalog_1$MAPPED_TRAIT_URI)
s_uni = lapply(s, function(x){
  uri = x$MAPPED_TRAIT_URI[1]
  print(uri)
  x$MAPPED_TRAIT = unname(efo$name[uri])
  return(x)
})
names(s_uni) = NULL
s_uni = do.call("c", s_uni)
catalog_1 = s_uni

# part 2: the trait that only contains letters
catalog_2 = GWAS_catalog[!(s1 - s2 > 0), ]

# make individual entries for snps linked with multiple phenotypes
catalog_2 = catalog_2 %>%
  rowwise() %>%
  mutate(MAPPED_TRAIT = str_split(MAPPED_TRAIT,pattern = ", "), MAPPED_TRAIT_URI = str_split(MAPPED_TRAIT_URI,pattern = ", ")) %>%
  unnest(c(MAPPED_TRAIT, MAPPED_TRAIT_URI)) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

catalog_2$MAPPED_TRAIT_URI = basename(catalog_2$MAPPED_TRAIT_URI)
catalog_2$MAPPED_TRAIT_URI = gsub("_", ":", catalog_2$MAPPED_TRAIT_URI)

catalog_all = c(catalog_2, catalog_1)
colnames(mcols(catalog_all))[c(35, 36)] = c("UNIQUE.TRAIT", "UNIQUE.URI") # to match the old name i have to work for follow up code

saveRDS(catalog_all, "./data/GWAS_catalog_noncoding_hg19_uniquetrait.RDS")
