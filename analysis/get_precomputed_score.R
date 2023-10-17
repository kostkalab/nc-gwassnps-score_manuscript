library(glmnet)
library(rtracklayer)
library(GenomicRanges)
library(stringr)

model_111diseases_DHS <- readRDS("./data/analysis/model_111diseases_DHS.RDS")

DIR_DHS = "./external_data/avocado/DHS/data/avocado.Dnase.bw/*.avocado.bw"
files = Sys.glob(DIR_DHS)
names = gsub(".DNase.allchrs.avocado.bw", "", basename(files))


chromosomes = paste0("chr", 1:22)
lapply(chromosomes, function(x){
  tmp = import(files[1], which=seqinfo(BigWigFile(files[1]))[x])
  smat = matrix(nrow = length(tmp), ncol = length(names))
  for (i in 1:127){
    tmp = import(files[i], which=seqinfo(BigWigFile(files[i]))[x])
    smat[,i] = tmp$score
    print(i)
  }
  colnames(smat) = names
  rm   = rowMeans(smat)
  X.n  = cbind(rm, smat)
  
  phenotypes = names(model_111diseases_DHS)
  preds = sapply(phenotypes, function(p){
    mdl = model_111diseases_DHS[[p]]
    pred = predict(mdl$glmnet.fit, 
                   newx = X.n, s = mdl$lambda.min, type = "response")
    towrite = tmp
    towrite$score = pred
    export.bw(towrite, str_c("./precomputed_score/",p,"_",x,".bw"))
    print(p)
    return(pred)
  })
  return(NULL)
})

