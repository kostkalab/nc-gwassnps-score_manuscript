
source("./model/fun_logreg_model.R")

get_prediction_score = function(PTPE, dat){
  message(PTPE)
  
  
  col         = colnames(dat)[grep("E[0-9]+", colnames(dat))]
  X           = as.matrix(dat[dat$phenotype==PTPE,col])
  rownames(X) = dat[dat$phenotype==PTPE,"snpID"]
  X[is.na(X)] = 0 #- FIXME
  y           = as.factor(dat$type[dat$phenotype==PTPE])
  
  rm   = rowMeans(X)
  X.c  = X # change this line
  X.n  = cbind(rm,X.c)
  
  # get prediction score
  set.seed(0525)
  mdl = mk_model(X=X.n, y=y)
  pred = predict(mdl$glmnet.fit, newx = X.n, s = mdl$lambda.min, type = "response")
  
  # take a look 
  y     = dat$type[dat$phenotype==PTPE]
  cat("Area under the PRC:", pr.curve(scores.class0=pred[,1],weights.class0=y)$auc.integral)
  print(roc(pred=pred[,1],resp=y,direction="<")$auc)
  
  # store prediction score
  pred_score = data.frame(snpID = rownames(X.n), 
                          phenotype = PTPE,
                          disease_associated = as.logical(y),
                          predicted_score = pred[,1])
  rownames(pred_score) = NULL
  
  return(list(mdl, pred_score))
}


