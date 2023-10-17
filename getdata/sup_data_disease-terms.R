library(readr)

# read in both coding and noncoding
GWAS_catalog_hg19 <- readRDS("./data/GWAS_catalog_hg19.RDS") # all SNVs

snpsnap_catalog_ct100_aggre <- readRDS("./data/catalog_cutoff100_snpsnap.RDS") 
catalog_ct100_aggre <- readRDS("./data/GWAS_catalog_filtered_cutoff100.RDS") # tell us whethet it is aggregated

phenos = sapply(catalog_ct100_aggre, function(x) x$AGGRE.TRAIT[1])
phenos = phenos[phenos %in% unique(snpsnap_catalog_ct100_aggre$phenotype)]



# get Number of SNVs used in the study
nsnps = sapply(phenos, function(x){
  d = snpsnap_catalog_ct100_aggre[snpsnap_catalog_ct100_aggre$type == 1]
  return(length(d[d$phenotype == x]))
})


# get Number of all SNVs in the GWAS catalog, including coding and noncoding SNVs
snps_all_nonaggre = sapply(names(phenos), function(x){
  nm = gsub(":", "_", x)
  b2 = GWAS_catalog_hg19[grep(nm, GWAS_catalog_hg19$MAPPED_TRAIT_URI),]
  b2 = unique(b2)
  return(length(b2))
})


# get Number of noncoding SNVs used in the study, before propagation
snps_noncoding_all_nonaggre = sapply(names(phenos), function(x){
  d1 = catalog_ct100_aggre[[x]]
  d1 = d1[!d1$aggre]
  d = snpsnap_catalog_ct100_aggre[snpsnap_catalog_ct100_aggre$type == 1]
  s = d[d$phenotype == d1$AGGRE.TRAIT[1]]
  return(length(d1[d1 %in% s]))
  
})


df_aggre = data.frame(term_name = phenos, 
                      term_efo_id = names(phenos), 
                      n_SNVs = snps_all_nonaggre, 
                      n_SNVs_noncoding_used = snps_noncoding_all_nonaggre, 
                      # n_SNVs_after_aggregation = nsnps_in_catalog, 
                      n_SNVS_noncoding_used_after_aggregation = nsnps)


df_aggre = df_aggre[order(df_aggre$term_name),]
rownames(df_aggre) = NULL

write_csv(df_aggre,"./sup_data/sup_data_disease-terms.csv.gz")


