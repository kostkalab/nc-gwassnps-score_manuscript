library(purrr)
library(dplyr)
library(GenomicRanges)

source("./getdata/fun_ld_marked.R")

## read in the file that contains clumped disease-associated SNVs
clumped_file = "./data/snpsnap_annotated/snpsnap_clumps_ld0.5/input_snps_clumped.txt"
clumped = read.table(clumped_file, sep = "\t", header = T, stringsAsFactors = F)

## read in the file for all GWAS SNVs, for only cutoff100 diseases. 
## Need this file because we need to extract p value for each SNV for each disease. 
####- Note: one SNV can have different p values in different diseases. so this should be done by each disease.
catalog_ct100_aggre <- readRDS("./data/GWAS_catalog_filtered_cutoff100.RDS")
names(catalog_ct100_aggre) = NULL
catalog_ct100_aggre = do.call("c", catalog_ct100_aggre) # unlist
chr = gsub("chr","",seqnames(catalog_ct100_aggre))
chr = gsub("X", "23", chr) # chromosome X -> chromosome 23
catalog_ct100_aggre$snpID = paste(chr, start(catalog_ct100_aggre), sep = ":") # annotate snpID

## read in the file for SNPSNAP matched SNVs, for cutoff100 diseases
snpsnap_catalog_ct100_aggre <- readRDS("./data/catalog_cutoff100_snpsnap.RDS")
snpsnap100 = snpsnap_catalog_ct100_aggre


ld_free_aggre = tapply(snpsnap100, snpsnap100$phenotype, ld_free)


# look at how many left
length(snpsnap100)
s1 = sort(tapply(snpsnap100, snpsnap100$phenotype, function(x){
  x = x[x$type == 1,]
  length(x)
})
)

s2 = sort(sapply(ld_free_aggre, function(x){
  x = x[x$type == 1 & x$ld_block_picked,]
  length(x)
})
)

df1 = data.frame(pheno = names(s1), original = s1)
df2 = data.frame(pheno = names(s2), ld_free = s2)
df = full_join(df1, df2, by = "pheno")
df$percent = df$ld_free/df$original

names(ld_free_aggre) = NULL
ld_free_aggre = do.call("c", ld_free_aggre)
saveRDS(ld_free_aggre, "./data/catalog_cutoff100_snpsnap_ld_marked.RDS")
