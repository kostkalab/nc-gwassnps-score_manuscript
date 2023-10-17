library(dplyr)
library(tibble)

coeff              <- readRDS("./data/analysis/coefficients.RDS")
top_mean           <- sapply(coeff, function(x) {x[,1]})
rownames(top_mean) <- rownames(coeff[[1]])
coeffi             <-  top_mean[-c(1,2),] # remove rowmeans and intercept, only concentrate on tissues
coeffi = t(coeffi)
coeffi = as_tibble(coeffi)
coeffi <- coeffi %>%
  add_column(phenotype = names(coeff),
             .before = "E001") 
readr::write_csv(coeffi, file="./sup_data/sup_data_beta-coefficients-mean-dhs.csv.gz") 


wts           <- sapply(coeff, function(x) {x[,2]})
rownames(wts) <- rownames(coeff[[1]])
wts           <- wts[-c(1,2),] # remove rowmeans and intercept, only concentrate on tissues
wts = t(wts)
wts = as_tibble(wts)
wts <- wts %>%
  add_column(phenotype = names(coeff),
             .before = "E001") 
readr::write_csv(wts, file="./sup_data/sup_data_beta-coefficients-sd-dhs.csv.gz") 
