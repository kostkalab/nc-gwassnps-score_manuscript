library(pROC)
library(ggplot2)
library(reshape2)
library(glmnet)
library(pamr)
library(ggbeeswarm)
library(PRROC)
library(BBmisc)
library(DT)
library(erma) # https://bioconductor.org/packages/devel/bioc/vignettes/erma/inst/doc/erma.html
library(gridExtra)
library(ggrepel)
library(uwot)
library(ComplexHeatmap)
library(circlize)
library(Rtsne)


source("./model/fun_logreg_get_coefficient.R")

# mk_model get from the source

# file to work with

DHS = readRDS("./data/anno_DHS_catalog_cutoff100_snpsnap.RDS")

phenotypes = tolower(gsub(" ","_",unique(DHS$phenotype)))


## get coefficients 
coeff = lapply(phenotypes, coeff_phenovs, dat = DHS, times = 10)
names(coeff) = phenotypes
saveRDS(coeff, "./data/analysis/coefficients.RDS")



