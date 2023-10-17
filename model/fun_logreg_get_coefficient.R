library(glmnet)



source("./model/fun_logreg_model.R")


# function to get the coefficients
coeff_phenovs = function(PTPE, dat, times = 10){
  # read in file. Each file contains a variant score to learn. 
  col         = colnames(dat)[grep("E[0-9]+", colnames(dat))]
  
  dat = dat[!is.na(rowSums(dat[,col])),]
  dat$phenotype = tolower(gsub(" ","_",dat$phenotype))
  
  # pick the disease to work with
  message(PTPE)
  
  # prepare X.n, y, rm for learning
  X           = as.matrix(dat[dat$phenotype==PTPE,col])
  X[is.na(X)] = 0 #- FIXME
  y           = as.factor(dat$type[dat$phenotype==PTPE])
  
  rm   = rowMeans(X)
  X.c  = X # change this
  X.n  = cbind(rm,X.c)
  
  # five fold. For each fold, learn training dataset and get coefficient. Five sets of coefficient
  get_coeff = function(X.n, y){
    flds = pamr:::balanced.folds(y,5)
    flds = lapply(flds,sort,decreasing=FALSE)
    
    # apply model on it and get predictions
    mkc <- function(fld){
      
      #- cv.glmnet result
      mdl  = mk_model(X=X.n[-fld,],y=y[-fld])
      # print(mdl$glmnet.fit$beta[1,1])
      print(mdl$lambda.min)
      coeff = as.matrix(coef(mdl, s = "lambda.min"))
      return(coeff)
    }
    
    coef_5 = lapply(flds,mkc)
    coef_5 = do.call("cbind", coef_5)
    return(coef_5)
  }
  
  # replicate 10 (default) times. 
  set.seed(0525)
  rep_coeff = replicate(times, get_coeff(X.n = X.n, y = y))
  m = rep_coeff[,,1] # make array into matrix
  for (i in 2:10){
    m = cbind(m, rep_coeff[,,i])
  }
  # now we have all 50 folds of value of coefficient for each tissue
  
  # get mean and sd for each tissue
  mn = apply(m,1,mean) # mean
  sd = apply(m,1,sd) # standard deviation
  df = data.frame(mean = mn, sd = sd)
  
  return(df)
}