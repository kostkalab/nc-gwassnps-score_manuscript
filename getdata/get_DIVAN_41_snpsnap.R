# get all matched snps
source("./getdata/fun_snpsnap_nmatch.R")

catalog_41           = readRDS("./data/divan_41.RDS")
catalog_all          = readRDS("./data/GWAS_catalog_noncoding_hg19_uniquetrait.RDS")
phegen_all           = readRDS("./data/phegen_processed_hg19.RDS")


# I already exclude all GWAS associated variants in the above function_snpsnap_nmatch.R file
g1KG_processed = g1KG_processed[!g1KG_processed %in% phegen_all] # exclude phegen variants in control variants
g1KG_processed = g1KG_processed[!g1KG_processed %in% catalog_all]


input_g       = "./data/snpsnap_annotated/divan_41_snpsnap/input_annotated.txt"
match_g       = "./data/snpsnap_annotated/divan_41_snpsnap/output_annotated.txt"

inputsnp_g    = read.table(input_g, header = T, stringsAsFactors = F, sep = "\t")
outputsnp_g   = read.table(match_g, header = T, stringsAsFactors = F, sep = "\t")


cat_matched = get_t_snpsnap_matched_snps(inputsnp_g, outputsnp_g)
cat_matched_df = data.frame(mcols(cat_matched), 
                            seqnames = seqnames(cat_matched), start = start(cat_matched), end = end(cat_matched))


# get matched snps i can use for catalog_41
catalog_41 = split(catalog_41, f = catalog_41$trait)
catalog_41_snpsnap = lapply(catalog_41, function(x){
  print(x$trait[1])
  g_coord = paste0(gsub("chr","",seqnames(x)),":", start(x))
  mc = data.frame(input_snp = g_coord, trait = x$trait, divan_new = x$divan_new)
  mc$input_snp = as.character(mc$input_snp)
  # move matched snps that are not in mc
  cat = cat_matched_df[cat_matched_df$input_snp %in% mc$input_snp, ]
  matched = full_join(mc, cat, by = "input_snp") # some snps don't have matched snps
  # remove snps that don't have matched snps
  matched = matched[!is.na(matched$snpID),]
  matched_snps = makeGRangesFromDataFrame(matched, keep.extra.columns = TRUE)
  names(matched_snps) = matched_snps$rsID
  return(matched_snps)
})
saveRDS(catalog_41_snpsnap, "./data/divan_41_snpsnap.RDS")

