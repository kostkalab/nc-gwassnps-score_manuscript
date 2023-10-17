library(GenomicRanges)
library(dplyr)
library(tibble)

catalog_41_snpsnap = readRDS("./data/divan_41_snpsnap.RDS")

# read in DIVAN score and DHS score
## read in DIVAN score
source("./getdata/fun_read_variantscore_organism.R")
source("./getdata/fun_read_variantscore_DIVAN.R")

catalog_DIVAN = lapply(catalog_41_snpsnap, function(x){
  dn = x$trait[1]
  gr = read_DIVAN(x, dn)
  return(gr)
})

## read in organism score
names(catalog_41_snpsnap) = NULL
catalog_41_snpsnap_unlist = do.call("c", catalog_41_snpsnap)
catalog_41_snpsnap_unlist = catalog_41_snpsnap_unlist[,c("trait", "divan_new", "type","snpID")] # only keep useful columns
catalog_org = read_GenoCanyon(catalog_41_snpsnap_unlist)

## read in DHS score
source("./getdata/fun_read_variantscore_tissue.R")
catalog_DHS = read_DHS(catalog_41_snpsnap_unlist)

names(catalog_DHS)[names(catalog_DHS) == "trait"] = "phenotype"


saveRDS(catalog_DIVAN, "./data/anno_divan41_DIVAN.RDS")
saveRDS(catalog_DHS, "./data/anno_divan41_DHS.RDS")
saveRDS(catalog_org, "./data/anno_divan41_GenoCanyon.RDS")