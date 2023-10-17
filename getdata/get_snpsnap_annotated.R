library(reshape2)


source("~/project/variant-scores_manuscript/getdata/fun_snpsnap_annotation.R")

output_file1 = "~/project/variant-scores_manuscript/snpsnap/v2_2020_04_27/output/GWAS_catalog_filtered_cutoff100_matched_snps_for_input_snps/matched_snps.tsv" # snpsnap
output_file2 = "~/project/variant-scores_manuscript/snpsnap/v2_2020_04_27/output/GWAS_catalog_filtered_cutoff100_tssonly_matched_snps_for_input_snps/matched_snps.tsv" # snpsnap-tss only

annotated_1 = snpsnap_annotation(output_file1)
annotated_2 = snpsnap_annotation(output_file2)

# check everything is good
summary(annotated_1[[1]]$dist_nearest_gene)
summary(annotated_1[[2]]$dist_nearest_gene)
summary(annotated_2[[2]]$dist_nearest_gene)

qqplot(x = annotated_1[[1]]$dist_nearest_gene, y = annotated_1[[2]]$dist_nearest_gene, main = "TSS qqplot: snpsnap 4 criteria")
qqplot(x = annotated_1[[1]]$dist_nearest_gene, y = annotated_2[[2]]$dist_nearest_gene, main = "TSS qqplot: snpsnap tss")


hist(annotated_1[[1]]$snp_maf)
hist(annotated_1[[2]]$snp_maf)
hist(annotated_2[[2]]$snp_maf)

summary(annotated_1[[1]]$gene_count)
summary(annotated_1[[2]]$gene_count)
summary(annotated_2[[2]]$gene_count)

summary(annotated_1[[1]]$friends_ld05)
summary(annotated_1[[2]]$friends_ld05)
summary(annotated_2[[2]]$friends_ld05)

write.table(annotated_1[[1]], file = "~/project/variant-scores_manuscript/data/snpsnap_annotated/catalog_filtered_cutoff100_snpsnap/input_annotated.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
write.table(annotated_1[[2]], file = "~/project/variant-scores_manuscript/data/snpsnap_annotated/catalog_filtered_cutoff100_snpsnap/output_annotated.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)

write.table(annotated_2[[1]], file = "~/project/variant-scores_manuscript/data/snpsnap_annotated/catalog_filtered_cutoff100_snpsnaptss/input_annotated.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
write.table(annotated_2[[2]], file = "~/project/variant-scores_manuscript/data/snpsnap_annotated/catalog_filtered_cutoff100_snpsnaptss/output_annotated.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
