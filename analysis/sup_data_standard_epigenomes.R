# retrieve ANATOMY, TYPE, and GROUP for 127 epigenomes. 
library(erma) # https://bioconductor.org/packages/devel/bioc/vignettes/erma/inst/doc/erma.html
library(readr)
meta          <- mapmeta()
met           <- as.data.frame(meta[,c("Epigenome.ID..EID.", "GROUP","Standardized.Epigenome.name", "ANATOMY", "TYPE")])
rownames(met) <- met[,1]
met = met[order(met$Epigenome.ID..EID.),]


readr::write_csv(met, file="./sup_data/sup_data_standard-epigenomes.csv.gz") 
