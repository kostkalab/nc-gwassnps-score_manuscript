library(GenomicRanges)
library(seqminer)
library(rtracklayer)

g1KG_processed              = readRDS("~/project/variant-scores_manuscript/data/1KG_processed_eur_noncoding.RDS") # snv, noncoding, eur

# add g_coord to g1kg.hg19 so i can match later on 
chr_g1KG = gsub("chr", "", seqnames(g1KG_processed))
g_cood.hg19 = paste(gsub("X","23",chr_g1KG),start(g1KG_processed), sep = ":")
g1KG_processed$g_cood.hg19 = g_cood.hg19


get_t_snpsnap_matched_snps = function(input, output, times = 10){
  a     = input
  b     = output
  # set.seed(0525) # to make this replicable
  set.seed(0123)
  
  ## clear up the data
  # random shuffle so duplicate deletion will be randomized between different input SNPs
  b     = b[sample(1:nrow(b)),] 
  # delete duplication
  b     = b[!duplicated(b$snpID),]
  
  # filter for only single nucleotide variants which is get from 1000 genome file that i previously generate
  b = b[b$snpID %in% g1KG_processed$g_cood.hg19,]
  
  # sample 10 matched snps from 100 matched snps, if less than 10, select them all
  
  t = tapply(b$snpID, b$input_snp, function(x) {
    if (length(x) > times){
      sa = sample(x, times)
    }else{
      sa = x
    }
  }) # 13128 out of 13135 input snps have 10 matched snps. 2 have 3 matched snps, 1 have 4 matched snps, ...
  # table(sapply(t, length))
  t = unlist(t)
  b_sampled = b[sort(match(t, b$snpID)),]
  
  ## make it grange object
  chr_a = sapply(strsplit(a$snpID,":"), function(x) x[[1]])
  pos_a = as.numeric(sapply(strsplit(a$snpID,":"), function(x) x[[2]]))
  chr_b = sapply(strsplit(b_sampled$snpID,":"), function(x) x[[1]])
  pos_b = as.numeric(sapply(strsplit(b_sampled$snpID,":"), function(x) x[[2]]))
  gr.a = GRanges(seqnames = chr_a, ranges = IRanges(start = pos_a, width = 1), type = rep(1, length(chr_a)))
  gr.b = GRanges(seqnames = chr_b, ranges = IRanges(start = pos_b, width = 1), type = rep(0, length(chr_b)))
  names(gr.a) = a$rsID
  names(gr.b) = b_sampled$rsID
  mcols(gr.a)[,2:4] = a[,c("snpID", "rsID", "snpID")]
  names(mcols(gr.a))[4] = c("input_snp")
  mcols(gr.b)[,2:4] = b_sampled[,c("snpID", "rsID", "input_snp")]
  gr = append(gr.a, gr.b)
  seqlevelsStyle(gr) <- "UCSC"
  return(gr)
}

