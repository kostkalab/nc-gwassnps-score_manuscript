library(GenomicRanges)
library(seqminer)
library(rtracklayer)



# 1 read in ARB and 1000 genome SNPs
catalog_aggre_snpsnap    = readRDS("./data/catalog_control_snpsnap.RDS")
g1KG_processed           = readRDS("./data/1KG_processed_eur_noncoding.RDS")
catalog_all              = readRDS("./data/GWAS_catalog_noncoding_hg19.RDS")
tss_19                   = readRDS("./data/tss_19.RDS")

catalog_aggre = catalog_aggre_snpsnap[catalog_aggre_snpsnap$type == 1]

# keep only autosome, that's what snpsnap did
g1KG_processed           = g1KG_processed[seqnames(g1KG_processed) %in% seqnames(catalog_aggre)] # keep the same chromosome as snpsnap did

# delete GWAS associated snps in the GWAS catalog from control snps
g1KG_processed           = g1KG_processed[!g1KG_processed %in% catalog_all] 

# get distance to nearest TSS
catlog.dist             = distanceToNearest(catalog_aggre, tss_19) # the catalog hg19
catalog_aggre$dist      = mcols(catlog.dist)[,1]


# 2 function to get TSS matched control SNPs
## cite function to get control snps using bins in control snps
source("./getdata/fun_get_control_snps_tss.R")

catalog_aggre_tss = get_control_snps_controlbins(catalog_aggre, g1KG_processed, bin = 50, times = 10) # all control snps
catalog_aggre_tss = c(catalog_aggre, catalog_aggre_tss) # disease and control snps
 
saveRDS(catalog_aggre_tss, "./data/catalog_control_tss.RDS")

