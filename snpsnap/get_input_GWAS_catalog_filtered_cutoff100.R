library(tibble)
library(stringr)
library(dplyr)
library(tidyr)

catalog_100_aggre_d = readRDS("~/project/variant-scores_manuscript/data/GWAS_catalog_filtered_cutoff100.RDS")

# generate a file that contains unique snps to match using snpsnap
cata_unlist = catalog_100_aggre_d
names(cata_unlist) = NULL
cata_unlist = do.call("c", cata_unlist)
chr = gsub("chr","",seqnames(cata_unlist))
chr = gsub("X", "23", chr)

g.coord = paste0(chr,":",start(cata_unlist))
g.coord = c("lead_snp",unique(g.coord))

# delete y chromosome
g.coord = g.coord[!grepl("Y", g.coord)]
write.table(g.coord, 
            file = "~/project/variant-scores_manuscript/snpsnap/input_GWAS_catalog_filtered_cutoff100.txt", 
            sep = "\n", row.names = F, col.names = F, quote = F)
