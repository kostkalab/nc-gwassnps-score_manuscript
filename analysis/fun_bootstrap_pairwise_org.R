library(sampler)
library(dplyr)
library(PRROC)
library(pROC)


bootstrap_samples = readRDS("./data/analysis/bootstrap_samples_30_0.9.RDS")



bootstrap_auc = function(phenotype, org, sc){
  print(phenotype)
  org_p = org[org$phenotype == phenotype]
  
  df = as.data.frame(mcols(org_p)[c(sc, "snpID", "type")])
  bootstrap_PTPE = data.frame(bootstrap_samples[[phenotype]])
  names(bootstrap_PTPE) = 1:ncol(bootstrap_PTPE)
  
  bootstrap_auc = lapply(bootstrap_PTPE, function(x){
    sub_df = df[df$snpID %in% x, ]
    # sub_df = sampler::ssamp(df=df, n=floor(nrow(df) * percent), strata=type) # sample 
    # sub_df = data.frame(sub_df)
    type = sub_df$type
    score = sub_df[,sc]
    
    ## for roc
    auroc = integer(length(score))
    for (i in seq_along(score)){
      auroc[i] = roc(type, score[,i])$auc
    }
    auc_roc = data.frame(PTPE = phenotype, auc = auroc, variant_score = sc)
    auc_roc$measure = "roc"
    
    ## for pr curve
    aupr = integer(length(score))
    for (i in seq_along(score)){
      s = !is.na(score[,i])
      aupr[i] = pr.curve(scores.class0=score[s,i],weights.class0=type[s])$auc.integral
    }
    auc_pr = data.frame(PTPE = phenotype, auc = aupr, variant_score = sc)
    auc_pr$measure = "pr"
    auc = rbind(auc_roc, auc_pr)
    return(auc)
  })
  for(i in seq_along(bootstrap_auc)){
    bootstrap_auc[[i]]$id = i
  }
  
  a = do.call("rbind", bootstrap_auc)
  return(a)
}
