library(GenomicRanges)

catalog_ct100_aggre     <- readRDS("./data/GWAS_catalog_filtered_cutoff100.RDS")
catalog_all_snpsnap     <- readRDS("./data/catalog_control_snpsnap.RDS")
catalog_all_snpsnap_tss <- readRDS("./data/catalog_control_snpsnap_tss.RDS")
catalog_all_tss         <- readRDS("./data/catalog_control_tss.RDS")
catalog_all_random      <- readRDS("./data/catalog_control_random.RDS")


## distributed matched snvs to each disease. 
get_match_SNPs = function(input_per_d, name, matched_all){
  # input_per_d tells us what disease snvs are in each disease
  # matched_all gives us matched control snvs for each disease snv
  chr = gsub("chr","",seqnames(input_per_d))
  chr = gsub("X", "23", chr)
  g_coord = paste0(chr, ":", start(input_per_d))
  input_per_d$snpsID = g_coord
  g = input_per_d[,c("SNPS", "snpsID", name)]
  mat = matched_all[matched_all$input_snp %in% g$snpsID]
  m = match(mat$input_snp, g$snpsID)
  mcols(mat)$phenotype = mcols(g)[m, c(name)]
  return(mat)
}


get_control_snps_for_each_disease_with_cutoff = function(disease_file, matched_file, name, cutoff = 100){
  catalog_cutoff100_bydisease = lapply(disease_file, get_match_SNPs, name = name, matched_all = matched_file)
  
  names(catalog_cutoff100_bydisease) = NULL
  
  catalog_cutoff100_bydisease = do.call("c", catalog_cutoff100_bydisease)
  
  ## filter again, only diseases with more than 100 snvs are kept. 
  a1 = catalog_cutoff100_bydisease[catalog_cutoff100_bydisease$type == 1]
  nms = names(table(a1$phenotype)[table(a1$phenotype) >= 100])
  catalog_cutoff100_bydisease = catalog_cutoff100_bydisease[catalog_cutoff100_bydisease$phenotype %in% nms]
  return(catalog_cutoff100_bydisease)
}


snpsnap_cutoff100     = get_control_snps_for_each_disease_with_cutoff(catalog_ct100_aggre, catalog_all_snpsnap,     name = "AGGRE.TRAIT", cutoff = 100)
snpsnaptss_cutoff100  = get_control_snps_for_each_disease_with_cutoff(catalog_ct100_aggre, catalog_all_snpsnap_tss, name = "AGGRE.TRAIT", cutoff = 100)
tss_cutoff100         = get_control_snps_for_each_disease_with_cutoff(catalog_ct100_aggre, catalog_all_tss,         name = "AGGRE.TRAIT", cutoff = 100)
random_cutoff100      = get_control_snps_for_each_disease_with_cutoff(catalog_ct100_aggre, catalog_all_random,      name = "AGGRE.TRAIT", cutoff = 100)

saveRDS(snpsnap_cutoff100,    "./data/catalog_cutoff100_snpsnap.RDS")
saveRDS(snpsnaptss_cutoff100, "./data/catalog_cutoff100_snpsnap_tss.RDS")
saveRDS(tss_cutoff100,        "./data/catalog_cutoff100_tss.RDS")
saveRDS(random_cutoff100,     "./data/catalog_cutoff100_random.RDS")

