library(pROC)
library(ggplot2)
library(reshape2)
library(glmnet)
library(pamr)
library(ggbeeswarm)
library(PRROC)

source("./model/fun_logreg_model.R")

train_and_test = function(X.train, y.train, X.test, y.test){
  mdl = mk_model(X=X.train,y=y.train)
  print(mdl$lambda.min)
  pred = predict(mdl$glmnet.fit,newx=X.test,s=c(mdl$lambda.min,mdl$lambda.1se,max(mdl$lambda)),type="response")
  
  ####
  # PR curve
  if (is.factor(y.test)) {
    y.test = as.integer(y.test)-1
  }
  # y = as.integer(y)-1
  r1 = pr.curve(scores.class0=pred[,1],weights.class0=y.test)$auc.integral
  r2 = pr.curve(scores.class0=pred[,2],weights.class0=y.test)$auc.integral
  r3 = pr.curve(scores.class0=pred[,3],weights.class0=y.test)$auc.integral
  res1 = c(r1,r2,r3)
  names(res1) = c("min","1se","max")
  res1 =c(res1, pr.curve(weights.class0=y.test,scores.class0=rowMeans(X.test))$auc.integral)
  names(res1)[4] = "baseline"
  
  #ROC curve
  r1 = roc(pred=pred[,1],resp=y.test,direction="<")$auc
  r2 = roc(pred=pred[,2],resp=y.test,direction="<")$auc
  r3 = roc(pred=pred[,3],resp=y.test,direction="<")$auc
  res2 = c(r1,r2,r3)
  names(res2) = c("min","1se","max")
  res2 =c(res2, roc(resp=y.test,pred=rowMeans(X.test), direction="<")$auc)
  names(res2)[4] = "baseline"
  
  ####
  res = cbind(res1, res2)
  colnames(res) = c("pr", "roc")
  return(res)
  
}


procPhenotype_by_chrom = function(PTPE, dat){
  #-------------------------------
  message(PTPE)
  
  col         = colnames(dat)[grep("E[0-9]+", colnames(dat))]
  traintest   = dat[dat$phenotype==PTPE,"chrom_heldout"]
  X           = as.matrix(dat[dat$phenotype==PTPE,col])
  X[is.na(X)] = 0 #- FIXME
  y           = as.factor(dat$type[dat$phenotype==PTPE])
  
  rm   = rowMeans(X)
  X.c  = X # change this line
  X.n  = cbind(rm,X.c)
  
  
  X.train = X.n[traintest == "train",]
  y.train = y[traintest == "train"]
  
  X.test = X.n[traintest == "test",]
  y.test = y[traintest == "test"]
  
  res = train_and_test(X.train, y.train, X.test, y.test)
  
  return(res)
}


procPhenotype_random = function(PTPE, dat){
  #-------------------------------
  message(PTPE)
  
  col         = colnames(dat)[grep("E[0-9]+", colnames(dat))]
  X           = as.matrix(dat[dat$phenotype==PTPE,col])
  X[is.na(X)] = 0 #- FIXME
  y           = as.factor(dat$type[dat$phenotype==PTPE])
  
  rm   = rowMeans(X)
  X.c  = X # change this line
  X.n  = cbind(rm,X.c)
  
  
  flds = pamr:::balanced.folds(y,5)
  # pick the first fold as randomly selected 1/5 test set. 
  ind = sort(flds[[1]])
  X.train = X.n[-ind,]
  y.train = y[-ind]
  
  X.test = X.n[ind,]
  y.test = y[ind]
  res2 = train_and_test(X.train, y.train, X.test, y.test)
  return(res2)
}


rep_procPheno_chr = function(PTPE, dat){
  set.seed(234324353)
  res_cd = procPhenotype_by_chrom(PTPE, dat = dat)
  df = melt(res_cd); colnames(df) = c("regularization","measure","auc")
  df$PTPE = PTPE
  return(df)
}
rep_procPheno_random = function(PTPE, times, dat){
  set.seed(234324353)
  res_cd = replicate(times,procPhenotype_random(PTPE, dat = dat))
  df = melt(res_cd); colnames(df) = c("regularization","measure","run_id","auc")
  df$PTPE = PTPE
  return(df)
}

