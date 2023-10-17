library(pROC)
library(ggplot2)
library(reshape2)
library(pamr)
library(ggbeeswarm)
library(PRROC)

source("./model/fun_logreg_model.R")

## cross-validation setup, and get the cross-validation auc
procPhenotype <- function(PTPE, dat){
  #-------------------------------
  message(PTPE)
  
  col         = colnames(dat)[grep("E[0-9]+", colnames(dat))]
  X           = as.matrix(dat[dat$phenotype==PTPE,col])
  X[is.na(X)] = 0 #- FIXME
  y           = as.factor(dat$type[dat$phenotype==PTPE])
  
  # set rowmeans
  rm   = rowMeans(X)
  X.c  = X # change this line
  X.n  = cbind(rm,X.c)
  
  flds = pamr:::balanced.folds(y,5) # five fold cross-validation
  flds = lapply(flds,sort,decreasing=FALSE)
  
  afu <- function(fld){
    #--------------------
    
    #- cv.glmnet result
    mdl  = mk_model(X=X.n[-fld,],y=y[-fld])
    # print(mdl$glmnet.fit$beta[1,1])
    print(mdl$lambda.min)
    pred = predict(mdl$glmnet.fit,newx=X.n[fld,],s=c(mdl$lambda.min,mdl$lambda.1se,max(mdl$lambda)),type="response")
    return(pred)
  }
  
  preds = lapply(flds,afu)
  
  ####
  
  ####
  
  y     = dat$type[dat$phenotype==PTPE]
  aucs  = sapply(1:5, function(x){
    r1 = pr.curve(scores.class0=preds[[x]][,1],weights.class0=y[flds[[x]]])$auc.integral
    r2 = pr.curve(scores.class0=preds[[x]][,2],weights.class0=y[flds[[x]]])$auc.integral
    r3 = pr.curve(scores.class0=preds[[x]][,3],weights.class0=y[flds[[x]]])$auc.integral
    res = c(r1,r2,r3)
    names(res) = c("min","1se","max")
    return(res)})
  res1 = rowMeans(aucs)
  res1 =c(res1, pr.curve(weights.class0=y,scores.class0=rowMeans(X))$auc.integral)
  names(res1)[4] = "baseline"
  
  #ROC curve
  aucs  = sapply(1:5, function(x){
    r1 = roc(pred=preds[[x]][,1],resp=y[flds[[x]]],direction="<")$auc
    r2 = roc(pred=preds[[x]][,2],resp=y[flds[[x]]],direction="<")$auc
    r3 = roc(pred=preds[[x]][,3],resp=y[flds[[x]]],direction="<")$auc
    res = c(r1,r2,r3)
    names(res) = c("min","1se","max")
    return(res)})
  
  res2 = rowMeans(aucs)
  res2 =c(res2, roc(resp=y,pred=rowMeans(X), direction="<")$auc)
  names(res2)[4] = "baseline"
  
  ####
  res = cbind(res1, res2)
  colnames(res) = c("pr", "roc")
  return(res)
}

# replicate 50 times
rep_procPheno = function(PTPE, times, dat){ 
  set.seed(234324353)
  res_cd = replicate(times,procPhenotype(PTPE, dat = dat))
  df = melt(res_cd); colnames(df) = c("regularization","measure","run_id","auc")
  df$PTPE = PTPE
  return(df)
}


