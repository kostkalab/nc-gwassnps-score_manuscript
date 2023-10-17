library(reshape2)


source("./getdata/fun_snpsnap_annotation.R")

output_file1 = "./snpsnap/v2_2020_04_27/output/divan_41_matched_snps_for_input_snps/matched_snps.tsv" # snpsnap

annotated_1 = snpsnap_annotation(output_file1)

# check everything is good
summary(annotated_1[[1]]$dist_nearest_gene)
summary(annotated_1[[2]]$dist_nearest_gene)

qqplot(x = annotated_1[[1]]$dist_nearest_gene, y = annotated_1[[2]]$dist_nearest_gene, main = "TSS qqplot: snpsnap 4 criteria")


hist(annotated_1[[1]]$snp_maf)
hist(annotated_1[[2]]$snp_maf)

summary(annotated_1[[1]]$gene_count)
summary(annotated_1[[2]]$gene_count)

summary(annotated_1[[1]]$friends_ld05)
summary(annotated_1[[2]]$friends_ld05)

write.table(annotated_1[[1]], file = "./data/snpsnap_annotated/divan_41_snpsnap/input_annotated.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
write.table(annotated_1[[2]], file = "./data/snpsnap_annotated/divan_41_snpsnap/output_annotated.txt", sep = "\t",
            row.names = F, col.names = T, quote = F)
