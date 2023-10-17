########
## this code is adapted from the following github repository: https://github.com/cboix/EPIMAP_ANALYSIS/blob/master/metadata_scripts/get_roadmap_mapping.R
#######


# Get list of files corresponding to 127 Roadmap Epigenomes:

library(jsonlite)
library(httr)
library(readxl)
today = format(Sys.time(), "%Y%m%d")

# -----------------------
# Get list of accessions:
# -----------------------
rdacc = c()
mainurl = 'https://www.encodeproject.org/'
urlpref = paste0(mainurl, 'report.tsv?type=Experiment&related_series.aliases=roadmap-epigenomics%3A')
# eidlist = scan('Annotation/roadmap_eidlist.tsv','c')

library(erma) # https://bioconductor.org/packages/devel/bioc/vignettes/erma/inst/doc/erma.html
meta          <- mapmeta()
met           <- as.data.frame(meta[,c("Epigenome.ID..EID.", "GROUP","Standardized.Epigenome.name", "ANATOMY", "TYPE")])
rownames(met) <- met[,1]

eidlist = met$Epigenome.ID..EID.


for (eid in eidlist){
  eidurl = paste0(urlpref, eid)
  call <- GET(eidurl) 
  df <- read.delim(text=rawToChar(call$content), skip=1)
  rdacc = rbind(rdacc, data.frame(eid = eid, Accession = as.character(df$Accession)))
}


data <- read_excel("./external_data/epimap/41586_2020_3145_MOESM5_ESM.xlsx", sheet = 2)
m = match(rdacc$Accession, data$Accession)

rdacc$bioiss = data$id[m]


unique_data <- rdacc %>%
  dplyr::select(eid, bioiss) %>%
  distinct() %>%
  filter(complete.cases(bioiss))

saveRDS(unique_data, "./external_data/epimap/roadmap_biosample_mapping.RDS")