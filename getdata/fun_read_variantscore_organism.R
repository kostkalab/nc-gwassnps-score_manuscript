library(seqminer)
library(rtracklayer)
library(GenomicRanges)
library(tibble)
library(dplyr)

# combine two granges, based on the first one

# make a tibble based on grange
# make_gr_tibble = function(g){
#   t1 = data.frame(seqnames=seqnames(g),
#                   start = start(g),
#                   end = end(g),
#                   snp= paste0(gsub("chr","",seqnames(g)), ":",start(g)), 
#                   mcols(g), stringsAsFactors = F)
#   return(t1)
# }
# combine_gr = function(gr, g_score){
#    t1 = make_gr_tibble(gr)
#    t2 = make_gr_tibble(g_score)[,-(1:3)]
#    df = full_join(t1,t2, by = "snp")
#    g = makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
# }


# readin CADD
read_CADD = function(gr.hg19, 
                     CADD_file = "./external_data/variant_scores/CADD/data/1000G_phase3.tsv.gz"){
  # check for duplicates, good. 01/27/20
  ARB = gr.hg19
  seqlevelsStyle(ARB) = "NCBI"
  ARB_Range = paste0(seqnames(ARB),":", start(ARB),"-", end(ARB))
  t = tabix.read.table(CADD_file, ARB_Range, col.names = TRUE,
                       stringsAsFactors = FALSE)
  t$cood.hg19 = paste0(t$Chrom, ":", t$Pos)
  t = t[!duplicated(t$cood.hg19),]
  CADD_g = GRanges(seqnames = t$Chrom, IRanges(start = t$Pos, width = 1))
  seqlevelsStyle(CADD_g) = "UCSC"
  f = findOverlaps(CADD_g, gr.hg19)
  g = granges(gr.hg19)
  mcols(g) = data.frame("RawScore" = NA, "PHRED" = NA)
  mcols(g[subjectHits(f)]) = t[queryHits(f),c("RawScore", "PHRED")]
  print("CADD")
  return(g)
}



# readin Eigen
# can read granges that contain duplicated snps, checked. 
read_eigen = function(gr.hg19, 
                      Eigen_dir = "./external_data/variant_scores/Eigen/data/",
                      type = "Eigen-PC"){
  Eigen_file = Sys.glob(paste(Eigen_dir,"*.tab.bgz",sep=""))
  ARB = gr.hg19
  seqlevelsStyle(ARB) = "NCBI"
  t_eigen = lapply(Eigen_file, function(x){
    chr = gsub("Eigen_hg19_noncoding_annot_chr|.tab.bgz", "", basename(x))
    gr = ARB[seqnames(ARB) == chr] 
    if (length(gr) > 0){
      gr_range = paste0(seqnames(gr),":", start(gr),"-", end(gr))
      t = tabix.read.table(x, gr_range, col.names = TRUE, 
                            stringsAsFactors = FALSE)
      t[,sapply(t,class) == "logical"] <- 
        sapply(t[,sapply(t,class) == "logical"], function(i) substr(as.character(i),1,1)) # in case some reads T as TRUE, change it back
      t2 = t[,c(1:4,31,33)]
      colnames(t2) = c("chr", "pos", "ref", "alt", "Eigen", "Eigen-PC")
      t2$id = paste0(t2$pos, ".", t2$alt)
      t2 = t2[!duplicated(t2$id),] # deal with duplicates
      t3 = t2[,c("chr","pos","alt",type)]
      t4 = t2[!duplicated(t2$pos),]
      t4 = t4[,c("chr","pos","ref",type)]
      colnames(t3)[3] = "allele"
      colnames(t4)[3] = "allele"
      t4[,type] = NA
      t5 = rbind(t3, t4)
      re = reshape(t5, v.names = type, idvar = "pos", timevar = "allele", direction = "wide")
      re$chr = paste0("chr", re$chr)
      l = c("A", "G", "T", "C")
      re = re[,c("chr", "pos", paste(type, l, sep = "."))]
      cat(x, "\n")
    } else{
      s = NULL
    }
    return(re)
  })
  eigen_score = do.call(rbind, t_eigen)
  eigen_score$mean = rowMeans(eigen_score[,3:6], na.rm = T)
  eigen_score_g = GRanges(seqnames = eigen_score$chr, IRanges(start = eigen_score$pos, width = 1))
  f = findOverlaps(eigen_score_g, gr.hg19)
  g = granges(gr.hg19)
  mcols(g) = data.frame("A" = NA, "G" = NA, "T" = NA, "C" = NA, mean = NA)
  mcols(g[subjectHits(f)]) = eigen_score[queryHits(f),3:7]
  return(g)
}


# read GenoCanyon
# works good if there are duplicated snps
read_GenoCanyon = function(gr.hg19, 
                           GenoCanyon_file = "./external_data/variant_scores/GenoCanyon/data/"){
  
  ARB = gr.hg19
  seqlevelsStyle(ARB) = "UCSC"
  
  IPF = Sys.glob(paste0(GenoCanyon_file,"*_hg19"))
  geno = lapply(IPF, function(x){
    files = Sys.glob(paste0(x,"/","*.tsv.gz"))
    chr = paste0("chr",gsub("GenoCanyon_Chr|_hg19", "", basename(x)))
    gr = ARB[seqnames(ARB) == chr]
    if(length(gr) > 0){
      gr_range = paste0(seqnames(gr),":", start(gr),"-", end(gr))
      sing = lapply(files, function(x){
        a = tabix.read.table(x, gr_range,  col.names = TRUE,
                             stringsAsFactors = FALSE)
        return(a)
      })
      sing = do.call("rbind", sing)
      colnames(sing) = c("chr", "pos", "score")
    } else{
      sing = NULL
    }
    return(sing)
  })
  geno = do.call("rbind", geno)
  geno = geno[,c(1,2,2,3)]
  names(geno)[2:4] = c("start", "end", "GenoCanyon")
  geno$snpID = paste0(gsub("chr", "", geno[,"chr"]), ":",geno[,"start"])
  g = geno[!duplicated(geno$snpID),]
  ARB = mcols(gr.hg19)
  
  a = full_join(as_tibble(ARB), g, by = "snpID")
  g = makeGRangesFromDataFrame(a, keep.extra.columns = TRUE)
  return(g)
}




# read GWAVA
# checked duplicated
read_GWAVA = function(gr.hg19, 
                      file = "./external_data/variant_scores/GWAVA/data/gwava_scores.bed.gz"){
  a = import.bed(file, which = gr.hg19, 
                 extraCols  = c(region = "numeric", TSS = "numeric", unmatched = "numeric"))
  m = match(gr.hg19, a)
  m = m[!is.na(m)]
  a = a[m]
  # meaning of three scores checked from GWAVA website
  f = findOverlaps(a, gr.hg19)
  g2 = granges(gr.hg19)
  mcols(g2) = data.frame("region" = NA, "TSS" = NA, "unmatched" = NA)
  mcols(g2[subjectHits(f)]) = mcols(a)[queryHits(f),c("region", "TSS", "unmatched")] # there was a warning before, be aware
  g3 = g2
  mcols(g3) = cbind(mcols(gr.hg19), mcols(g2))
  return(g3)
}


# read linsight
# duplicated snps check 
read_LINSIGHT = function(gr.hg19, 
                         file = "./external_data/variant_scores/linsight/data/LINSIGHT.bw"){
  a = import(file, which = disjoin(gr.hg19))
  f = findOverlaps(a, gr.hg19)
  g = gr.hg19
  g$score = NA
  g$score[subjectHits(f)] = a[queryHits(f)]$score
  return(g)
}



# read ncER
read_ncER = function(gr.hg19, 
                     folder = "./external_data/variant_scores/ncER/data/"){
  IPF = Sys.glob(paste0(folder,"*coordSorted.txt.gz"))
  
  g = gr.hg19
  g$score = NA_real_
  
  for (i in seq_along(IPF)){
    print(IPF[i])
    gz = read.table(gzfile(IPF[i]),sep="\t")
    colnames(gz) = c("chromosome", "start", "end", "score")
    gz$start = gz$start + 1
    gr = makeGRangesFromDataFrame(gz
                                  , keep.extra.columns=T)
    f = findOverlaps(gr, gr.hg19)
    g$score[subjectHits(f)] = gr[queryHits(f)]$score
  }
    return(g)

}


read_CADD_wholegenome = function(gr.hg19, 
                                 CADD_file = "./external_data/variant_scores/CADD/data/v1.6_hg19/whole_genome_SNVs.tsv.gz"){
  # check for duplicates, good. 01/27/20
  ARB = gr.hg19
  seqlevelsStyle(ARB) = "NCBI"
  ARB_Range = paste0(seqnames(ARB),":", start(ARB),"-", end(ARB))
  t = tabix.read.table(CADD_file, ARB_Range, col.names = TRUE,
                       stringsAsFactors = FALSE)
  t$cood.hg19 = paste0(t$Chrom, ":", t$Pos)
  # get the mean score of three alleles
  mean_score = enframe(tapply(t$PHRED, t$cood.hg19, mean))
  colnames(mean_score) = c("cood.hg19", "mean_PHRED")
  # mean_score$chrom = gsub("([0-9]|[A-Z]):([0-9]+)", "\\1", mean_score$cood.hg19)
  # mean_score$pos  = gsub("([0-9]|[A-Z]):([0-9]+)", "\\2", mean_score$cood.hg19)
  t1 = full_join(t,mean_score, by = "cood.hg19")
  t1 = t1[!duplicated(t1$cood.hg19),]
  CADD_g = GRanges(seqnames = t1$Chrom, IRanges(start = t1$Pos, width = 1))
  seqlevelsStyle(CADD_g) = "UCSC"
  f = findOverlaps(CADD_g, gr.hg19)
  g = granges(gr.hg19)
  mcols(g) = data.frame("RawScore" = NA, "PHRED" = NA)
  mcols(g[subjectHits(f)]) = t1[queryHits(f),c("RawScore", "mean_PHRED")]
  return(g)
}
