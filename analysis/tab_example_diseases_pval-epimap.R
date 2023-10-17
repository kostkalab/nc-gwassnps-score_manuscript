library(dplyr)
library(xtable)
library(readr)
library(readxl)

# Run "./external_data/epimap/get_roadmap_biosample_mapping.R" to get "./external_data/epimap/roadmap_biosample_mapping.RDS" file
mapping = readRDS("./external_data/epimap/roadmap_biosample_mapping.RDS")
gwas_results_hg_snpcentric = readRDS("./external_data/epimap/gwas_results_hg_snpcentric.Rds")
top5 = readRDS("./data/analysis/example_disease_top5_tissues.RDS")

format_p_values <- function(p_values) {
  formatted_p_values <- ifelse(p_values < 0.001, "<0.001",
                               ifelse(p_values <= 0.01, round(p_values, digits = 3),
                                      round(p_values, digits = 2)))
  return(formatted_p_values)
}



a1 = gwas_results_hg_snpcentric[grep("Systemic sclerosis" ,gwas_results_hg_snpcentric$uid,  ignore.case = TRUE),]
table(a1$uid)
a1_sig = a1 %>% filter(padj < 0.1)


m = match( a1_sig$id, mapping$bioiss)
a1_sig$eid = mapping$eid[m]


p1 = sapply(top5$systemic_scleroderma$Epigenome.ID..EID., function(eid){
  m = which(!is.na(match(mapping$eid, eid)) )
  bioss = mapping$bioiss[m]
  padj = a1[a1$id %in% bioss,]
  p = format_p_values(padj$padj[1])
})

p1

t1 = top5$systemic_scleroderma
# t1$padj =

a2 = gwas_results_hg_snpcentric[grep("cholangitis" ,gwas_results_hg_snpcentric$uid,  ignore.case = TRUE),]
table(a2$uid)

p2 = sapply(top5$sclerosing_cholangitis$Epigenome.ID..EID., function(eid){
  m = which(!is.na(match(mapping$eid, eid)) )
  bioss = mapping$bioiss[m]
  padj = a2[a2$id %in% bioss,]
  p = format_p_values(padj$padj[1])
})

p2
a3 = gwas_results_hg_snpcentric[grep("Colorectal" ,gwas_results_hg_snpcentric$uid, ignore.case = TRUE),]
table(a3$uid)

p3 = sapply(top5$colorectal_adenoma$Epigenome.ID..EID., function(eid){
  m = which(!is.na(match(mapping$eid, eid)) )
  bioss = mapping$bioiss[m]
  padj = a3[a3$id %in% bioss,]
  p = format_p_values(padj$padj[1])
})

p3

a4 = gwas_results_hg_snpcentric[grep("Atrial fibrillation" ,gwas_results_hg_snpcentric$uid, ignore.case = TRUE),]
table(a4$uid)

p4 = sapply(top5$atrial_fibrillation$Epigenome.ID..EID., function(eid){
  m = which(!is.na(match(mapping$eid, eid)) )
  bioss = mapping$bioiss[m]
  padj = a4[a4$id %in% bioss,]
  p = format_p_values(padj$padj[1])
})

p4
a5 = gwas_results_hg_snpcentric[grep(" Cutaneous " ,gwas_results_hg_snpcentric$uid,  ignore.case = TRUE),]
table(a5$uid)

p5 = sapply(top5$cutaneous_melanoma$Epigenome.ID..EID., function(eid){
  m = which(!is.na(match(mapping$eid, eid)) )
  bioss = mapping$bioiss[m]
  padj = a5[a5$id %in% bioss,]
  p = format_p_values(padj$padj[1])
})
p5


p1
p2
p3
p4
p5

