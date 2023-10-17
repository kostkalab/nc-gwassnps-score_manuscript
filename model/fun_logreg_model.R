library(glmnet)

#- fits / assesses l2-regulrized LR model
mk_model <- function(X,y){
  #=========================
  
  #-observation weights
  w       = rep(1/sum(y==0),length(y))
  w[y==1] = 1/sum(y==1)
  w       = w/max(w)
  
  #- exclude first feature from penalty
  pf = c(0,rep(1,ncol(X)-1))
  
  res = cv.glmnet(X, y, nfolds=5,  type.measure="auc", family="binomial", alpha=0, weights=w, intercept=F, standardize=FALSE,penalty.factor=pf, lower.limits = c(0, rep(-Inf, ncol(X-1))))
  
  return(res)
}
