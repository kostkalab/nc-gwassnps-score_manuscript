library(GenomicRanges)
library(seqminer)
library(rtracklayer)
library(ggplot2)


# 1 read in disease-associated and 1000 genome SNPs
catalog_aggre_snpsnap    = readRDS("./data/catalog_control_snpsnap.RDS")
g1KG_processed           = readRDS("./data/1KG_processed_eur_noncoding.RDS")
catalog_all              = readRDS("./data/GWAS_catalog_noncoding_hg19.RDS")

catalog_aggre = catalog_aggre_snpsnap[catalog_aggre_snpsnap$type == 1]

# keep only autosome, that's what snpsnap did
g1KG_processed           = g1KG_processed[seqnames(g1KG_processed) %in% seqnames(catalog_aggre)] # keep the same chromosome as snpsnap did

# delete gwas catalog from g1KG processed, double check before use. Add this line on June 2, 2020, after i changed g1kg_processed to include ARB SNPs
g1KG_processed           = g1KG_processed[!g1KG_processed %in% catalog_all] 


# 3 function to get random controls
set.seed(0525)
r_matched = sample(g1KG_processed, length(catalog_aggre)*10) # all control snps
r_matched2             = GRanges(granges(r_matched), 
                            type = 0, 
                            snpID = paste(gsub("chr","",seqnames(r_matched)), start(r_matched), sep = ":"), 
                            rsID = names(r_matched), 
                            input_snp = rep(catalog_aggre$input_snp, each = 10)) # assigned inputsnps to them - it's already random

catalog_aggre_random = c(catalog_aggre, r_matched2) # disease and control snps

saveRDS(catalog_aggre_random, "./data/catalog_control_random.RDS")


