require(lpSolve)
library(dplyr)


get_set <- function(vpc,frac){
  #=============================
  
  n_chroms <- nrow(vpc)
  w_pos    <- vpc[,"pos"] / sum(vpc[,"pos"])
  w_neg    <- vpc[,"neg"] / sum(vpc[,"neg"])
  
  f.obj <- w_pos - abs(w_pos - w_neg)
  f.con <- t(w_pos)
  f.dir <- c("<=")
  f.rhs <- c(frac)
  
  sol  <- lp("max", f.obj, f.con, f.dir, f.rhs, binary.vec = 1:n_chroms)
  inds <- which(sol$solution ==1)
  
  rr <- list(inds     = inds,
             frac_pos = sum(w_pos[inds]),
             frac_neg = sum(w_neg[inds]))
  return(rr)
}

get_chrom_heldout = function(phenotype, matched_file){
  print(phenotype)
  gr_pheno = matched_file[matched_file$phenotype == phenotype] # get the variants for that disease
  gr_pheno$chrom = factor(gsub("chr","",seqnames(gr_pheno)), levels = 1:22)
  gr1 = gr_pheno[gr_pheno$type == 1] # get all disease associated variants
  gr2 = gr_pheno[gr_pheno$type == 0] # get all control variants
  
  df1 = data.frame(table(gr1$chrom))
  names(df1) = c("chr", "pos")
  df2 = data.frame(table(gr2$chrom))
  names(df2) = c("chr", "neg")
  df3 = full_join(df1, df2, by = "chr")
  
  hold_out <- get_set(df3, frac = 0.2) # hold out ind is the chromosome
  gr_pheno$chrom_heldout = "train"
  gr_pheno[gr_pheno$chrom %in% hold_out$inds]$chrom_heldout = "test"
  return(gr_pheno)
}
