library(purrr)
library(dplyr)
library(GenomicRanges)


give_min = function(a){
  b = integer(length(a))
  b[which.min(a)] = 1
  return(b)
}

## this function is the main function in this code file. 
## input: 1.file to be matched: catalog_cutoff100_xxx.RDS; 
## input: 2.clumped file to give us ld cluster; 
## input: 3.gwas catalog file to give us p.value
ld_free = function(x, clump_file = clumped, gwas_file = catalog_ct100_aggre){
  ### for each disease
  print(x$phenotype[1])
  inputsnp = unique(x$input_snp)
  
  clump_file$clump_file_snps = paste(clump_file$clumped_snps, ";", sep = "")
  clump_file$combined = paste(clump_file$index_snp, clump_file$clumped_snps, sep = ";") # i add ; to the end in case i search for 1:12345 and there is 1:123456, which is not it
  s = sapply(inputsnp, function(x){
    snp = paste0(x, ";")
    grep(snp, clump_file$combined)
    # give each SNV a cluster id. The cluster id is the rank of this snv in the clumped file. 
    # So same SNV will have the same cluster id, thus we cluster the SNVs into LD blocks.  
  })
  
  cd = gwas_file[gwas_file$AGGRE.TRAIT == x$phenotype[1]] # get all catalog cutoff100 variants for this disease
  df = data.frame(snpID = names(s), cluster = s)
  df$snpID = as.character(df$snpID)
  rownames(df) = NULL
  df_p = as.data.frame(mcols(cd)[,c("P.VALUE", "snpID")]) # get the p values for all disease-associated SNVs
  df = left_join(df, df_p, by = "snpID") # join together by snpID. 
  tib <- tibble(df) %>%
    group_by(cluster) %>%
    mutate(ld_block_picked = give_min(P.VALUE)) # find the one with the smallest p.value, mark it as ld_block_picked 
  
  # map back
  m = match(x$input_snp, tib$snpID)
  gr = x
  mcols(gr) = cbind(mcols(gr), tib[m,2:4])
  return(gr)
}

