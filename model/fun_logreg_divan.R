library(pROC)
library(ggplot2)
library(reshape2)
library(glmnet)
library(pamr)
library(ggbeeswarm)
library(PRROC)

test_boostrap = readRDS("./data/analysis/bootstrap_samples_divan29.RDS")
pheno = tolower(gsub(" ","_",names(test_boostrap)))
names(test_boostrap) = pheno

source("./model/fun_logreg_model.R")
#  #- fits / assesses l2-regulrized LR model
#  mk_model <- function(X,y){
#    #=========================
#    
#    #-observation weights
#    w       = rep(1/sum(y==0),length(y))
#    w[y==1] = 1/sum(y==1)
#    w       = w/max(w)
#    
#    #- exclude first feature from penalty
#    pf = c(0,rep(1,ncol(X)-1))
#    
#    #- assess lambdas via 5-fold CV
#    # res = cv.glmnet(X, y, nfolds=5, family="binomial", alpha=0, weights=w, intercept=TRUE, standardize=FALSE,penalty.factor=pf)
#    # res = cv.glmnet(X, y, nfolds=5,  type.measure="auc", family="binomial", alpha=0, lambda = 2^seq(4,-12, length.out = 100), weights=w, intercept=TRUE, standardize=FALSE,penalty.factor=pf)
#    res = cv.glmnet(X, y, nfolds=5,  type.measure="auc", family="binomial", alpha=0, weights=w, intercept=F, standardize=FALSE,penalty.factor=pf, lower.limits = 0)
#    return(res)
#  }
#  

procPhenotype_oldnew <- function(PTPE, dat){
  #-------------------------------
  set.seed(0525)
  bts = data.frame(test_boostrap[[PTPE]], stringsAsFactors = F)
  
  message(PTPE)
  
  X         = as.matrix(dat[dat$phenotype==PTPE,seq_len(ncol(dat)-4)])
  X[is.na(X)] = 0 #- FIXME
  y         = as.factor(dat$type[dat$phenotype==PTPE])
  divan_new = dat$divan_new[dat$phenotype == PTPE]
  snpID     = dat$snpID[dat$phenotype == PTPE]
  
  
  rm   = rowMeans(X)
  X.c  = X - rm
  X.n  = cbind(rm,X.c)
  
  # seperate training and testing dataset
  # training dataset
  # X_training
  # y_training
  
  X_training = X.n[divan_new == 0,]
  y_training = y[divan_new == 0]
  
  
  mdl = mk_model(X=X_training, y=y_training)
  # print(mdl$glmnet.fit$beta[1,1])
  print(mdl$lambda.min)
  
  
  # testing dataset
  # X_test
  X_test = X.n[divan_new == 1,]
  # y_test
  y_test = y[divan_new == 1]
  id_test = snpID[divan_new == 1]
  
  # test make prediction
  pred = predict(mdl$glmnet.fit,newx=X_test,s=c(mdl$lambda.min,mdl$lambda.1se,max(mdl$lambda)),type="response")
  pred = as.data.frame(pred)
  pred$snpID = id_test
  names(pred)[1:3] = c("min", "1se", "max")
  pred$type = as.numeric(y_test) -1
  
  # measure test dataset performance
  
  auc1 = sapply(bts, function(x){
    
    X_samp = pred[pred$snpID %in% x,]
    auroc = roc(X_samp$type, X_samp$min)$auc
    aupr = pr.curve(scores.class0=X_samp$min,weights.class0=X_samp$type)$auc.integral
    res1 = c(auroc, aupr)
    names(res1) = c("roc", "pr")
    return(res1)
  })
  aucd = reshape2::melt(auc1)
  colnames(aucd) = c("type", "id", "auc")
  aucd$PTPE = PTPE
  
  return(aucd)
}
