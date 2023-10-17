library(GenomicRanges)
library(seqminer)
library(rtracklayer)
library(ggplot2)
library(tibble)
library(dplyr)
source("./getdata/fun_read_variantscore_organism.R")

# notice: there should be no duplication in grange object when read in variant scores
# remove this notice


# read in file

catalog_snpsnap_aggre       <- readRDS("./data/catalog_cutoff100_snpsnap.RDS")
catalog_snpsnap_tss_aggre    <- readRDS("./data/catalog_cutoff100_snpsnap_tss.RDS")
catalog_tss_aggre           <- readRDS("./data/catalog_cutoff100_tss.RDS")
catalog_random_aggre        <- readRDS("./data/catalog_cutoff100_random.RDS")


# store out even if one of them is not working. 
read_organism_score = function(g){
  g1 = g 
  # g1 = split(g1, f = g1$phenotype)
  g1$CADD       = read_CADD  (g)$PHRED
  # eigenscore = tapply(g1, g1$phenotype, read_eigen)
  g1$eigen      = read_eigen (g)$mean
  g1$GenoCanyon = read_GenoCanyon(g)$GenoCanyon
  g1$GWAVA      = read_GWAVA(g)$TSS
  g1$LINSIGHT   = read_LINSIGHT(g)$score
  return(g1)
}


catalog_snpsnap_aggre_org         = read_organism_score(catalog_snpsnap_aggre)
catalog_tss_aggre_org             = read_organism_score(catalog_tss_aggre)

saveRDS(catalog_snpsnap_aggre_org ,     "./data/anno_org_catalog_cutoff100_snpsnap.RDS")
saveRDS(catalog_tss_aggre_org ,         "./data/anno_org_catalog_cutoff100_tss.RDS")



catalog_snpsnap_tss_aggre_org     = read_organism_score(catalog_snpsnap_tss_aggre)
catalog_random_aggre_org          = read_organism_score(catalog_random_aggre)

saveRDS(catalog_snpsnap_tss_aggre_org,  "./data/anno_org_catalog_cutoff100_snpsnaptss.RDS")
saveRDS(catalog_random_aggre_org,       "./data/anno_org_catalog_cutoff100_random.RDS")