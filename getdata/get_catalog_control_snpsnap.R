library(GenomicRanges)
library(seqminer)
library(rtracklayer)


source("~/project/variant-scores_manuscript/getdata/[fun]snpsnap_nmatch.R")
GWAS_catalog_noncoding_hg19_uniquetrait <- readRDS("~/project/variant-scores_manuscript/data/GWAS_catalog_noncoding_hg19_uniquetrait.RDS")
g1KG_processed = g1KG_processed[!g1KG_processed %in% GWAS_catalog_noncoding_hg19_uniquetrait] # exclude GWAS associated variants from control snps

inputsnp1 = read.table("~/project/variant-scores_manuscript/data/snpsnap_annotated/catalog_filtered_cutoff100_snpsnap/input_annotated.txt", 
                       header = T, stringsAsFactors = F, sep = "\t")
outputsnp1 = read.table("~/project/variant-scores_manuscript/data/snpsnap_annotated/catalog_filtered_cutoff100_snpsnap/output_annotated.txt", 
                        header = T, stringsAsFactors = F, sep = "\t")
inputsnp2 = read.table("~/project/variant-scores_manuscript/data/snpsnap_annotated/catalog_filtered_cutoff100_snpsnaptss/input_annotated.txt", 
                       header = T, stringsAsFactors = F, sep = "\t")
outputsnp2 = read.table("~/project/variant-scores_manuscript/data/snpsnap_annotated/catalog_filtered_cutoff100_snpsnaptss/output_annotated.txt", 
                        header = T, stringsAsFactors = F, sep = "\t")

gr_cat1 = get_t_snpsnap_matched_snps(inputsnp1, outputsnp1)
gr_cat2 = get_t_snpsnap_matched_snps(inputsnp2, outputsnp2)

# there are some snps not matched by snpsnap but are matched by snpsnap tss
a1 = gr_cat2[gr_cat2$type == 1] # 25851
a0 = gr_cat2[gr_cat2$type == 0]

table(a1 %in% a0) # check if there's any input snps in output snps
a1 = a1[a1 %in% gr_cat1] # 25540
a0 = a0[a0$input_snp %in% a1$input_snp] # 221520

gr_cat2 = c(a1, a0) # keep the same input variants for snpsnap for snpsnap_tss so they are comparable

# s = split(gr_cat, f = gr_cat$input_snp)
# table(sapply(s, length)) # 22145 out of 22157 have matched with 10 snps

# saveRDS(gr_cat, "~/projects/variant-scores-2/data/GWAS_catalog_project_3.0/genetic/catalog_all_snpsnap.RDS")

saveRDS(gr_cat1, "~/project/variant-scores_manuscript/data/catalog_control_snpsnap.RDS")
saveRDS(gr_cat2, "~/project/variant-scores_manuscript/data/catalog_control_snpsnap_tss.RDS")