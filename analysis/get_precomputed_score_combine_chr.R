library(rtracklayer)
library(GenomicRanges)
library(stringr)

# combine chromosomes
files = Sys.glob("~/rstudio_srv/project/variant-scores_manuscript/precomputed_score/*")
files_base = basename(files)
phenotypes = unique(gsub("_chr[0-9]*.bw", "", files_base))
lapply(phenotypes, function(p){
  print(p)
  file = Sys.glob(paste0("~/rstudio_srv/project/variant-scores_manuscript/precomputed_score/",p,"_*"))
  bws = lapply(file, import)
  bws = do.call("c", bws)
  export.bw(bws, str_c("~/rstudio_srv/project/variant-scores_manuscript/precomputed_score_combined/",p,".bw"))
})