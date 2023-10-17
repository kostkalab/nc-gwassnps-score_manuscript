library(dplyr)

# get coefficients and epigenome annotation
#------------------------------------------

coeff              <- readRDS("./data/analysis/coefficients.RDS")
#- get mean of coefficients
top_mean           <- sapply(coeff, function(x) {x[,1]})
rownames(top_mean) <- rownames(coeff[[1]])
coeffi             <-  top_mean[-c(1,2),] # remove rowmeans and intercept, only concentrate on tissues


# retrieve ANATOMY, TYPE, and GROUP for 127 epigenomes. 
library(erma) # https://bioconductor.org/packages/devel/bioc/vignettes/erma/inst/doc/erma.html
meta          <- mapmeta()
met           <- as.data.frame(meta[,c("Epigenome.ID..EID.", "GROUP","Standardized.Epigenome.name", "ANATOMY", "TYPE")])
rownames(met) <- met[,1]
met           <- met[rownames(coeffi),]



#- weights are bootstrap variances
wts           <- sapply(coeff, function(x) {x[,2]})
rownames(wts) <- rownames(coeff[[1]])
wts           <- wts[-c(1,2),] # remove rowmeans and intercept, only concentrate on tissues

get_wghtd_corrs <- function(X,wts){
  md = quantile(wts,.25)/4
  get_wghtd_cor <- function(i,j,omit,signed=TRUE){
    twts =  1/sqrt(wts[     ,i]*3/4+md)*1/sqrt(wts[     ,j]*3/4+md)
    tmp <- lm(X[     ,i] ~ X[     ,j], weights =twts)
    if(!signed) return(summary(tmp)$r.squared)
    return(summary(tmp)$r.squared * sign(tmp$coefficients[2]))
  }
  res = matrix(NA,ncol=ncol(X),nrow=ncol(X))
  for(i in 2:ncol(X)){
    for(j in 1:(i-1)){
      res[i,j] = get_wghtd_cor(i,j,omit=omit)
      res[j,i] = res[i,j]
    }
  }
  diag(res) = 1
  colnames(res) = colnames(X)
  rownames(res) = colnames(X)
  return(res)
}

### get weighted correlations
tmp <- get_wghtd_corrs(coeffi,wts) # get weighted correlations
plot(cor(coeffi),tmp, pch=".") # compared correlations with weighted correlations


saveRDS(tmp, "./data/analysis/coeffi_weighted_correlations.RDS")

