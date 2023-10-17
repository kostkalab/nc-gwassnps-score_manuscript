# read DIVAN
# not sure if it can deal with duplicated snps
read_DIVAN = function(gr.hg19, diseasename){
  # setwd("~/project/variant-scores_manuscript/external_data/variant_scores/DIVAN/allTSSscore/")
  d = data.frame(seqnames(gr.hg19), start(gr.hg19), end(gr.hg19))
  write.table(d, file = "./external_data/variant_scores/DIVAN/tmpDisease.txt", sep = "\t", row.names = F, col.names = F, quote = F)
  queryfile="./external_data/variant_scores/DIVAN/tmpDisease.txt"
  disease=gsub(" ", "", diseasename, fixed = TRUE)
  disease = paste0("./external_data/variant_scores/DIVAN/allTSSscore/", disease)
  distribution="./external_data/variant_scores/DIVAN/allTSSscore/scoredistTSS"
  scorefile="./external_data/variant_scores/DIVAN/tmp-score.txt"
  source("./external_data/variant_scores/DIVAN/allTSSscore/DIVANtoolkit/scoreDIVAN.console.genome.R")
  scoreDIVAN(queryfile,disease,distribution,scorefile)
  DIVANscore = read.table(scorefile, header = T)
  DIVANscore$snpID = paste0(gsub("chr","",DIVANscore$chr), ":", DIVANscore$start)
  gr_df = data.frame(mcols(gr.hg19), 
                     seqnames = seqnames(gr.hg19), start = start(gr.hg19), end = end(gr.hg19))
  gr_divan = full_join(gr_df, DIVANscore[,-(1:3)], by = "snpID") # full join keep the one that contains more
  # gr_df should contain more snps than DIVANscore
  DIVAN_g = makeGRangesFromDataFrame(gr_divan, keep.extra.columns = TRUE)
  
  return(DIVAN_g)
}
