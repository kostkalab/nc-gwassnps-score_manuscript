library(dplyr)

DHS_SNPSNAP_logreg <- readRDS("./data/analysis/perf_DHS_logreg.RDS")

DHS_SNPSNAP_logreg %>% filter(regularization == "min" & measure == "pr") %>% 
  group_by(PTPE) %>%
  summarise(mean_auc = mean(auc)) %>% 
  arrange(desc(mean_auc)) %>%
  print(n = 111)





coeff              <- readRDS("./data/analysis/coefficients.RDS")
#- get mean of coefficients
top_mean           <- sapply(coeff, function(x) {x[,1]})
rownames(top_mean) <- rownames(coeff[[1]])
coeffi             <-  top_mean[-c(1,2),] # remove rowmeans and intercept, only concentrate on tissues


library(erma) # https://bioconductor.org/packages/devel/bioc/vignettes/erma/inst/doc/erma.html
meta          <- mapmeta()
met           <- as.data.frame(meta[,c("Epigenome.ID..EID.", "GROUP","Standardized.Epigenome.name", "ANATOMY", "TYPE")])
rownames(met) <- met[,1]
met           <- met[rownames(coeffi),]


examples = c("systemic_scleroderma", "sclerosing_cholangitis", "colorectal_adenoma", 
             # "neoplasm_of_mature_b-cells", 
             #"adult_onset_asthma", 
             "atrial_fibrillation",
             "cutaneous_melanoma"
             # , 
             # "hypothyroidism", 
             #"alzheimer's_disease"
             )
coeffi_examples = coeffi[,examples]




top5 = apply(coeffi_examples, 2, function(x){
  nms = names(x)
  top10 = names(sort(x, decreasing = T)[1:5])
  nms = met[top10,]
  nms_tab = nms[,c("Epigenome.ID..EID.","Standardized.Epigenome.name", "ANATOMY")]
  nms_tab$rank = 1:5
  nms_tab = nms_tab[,c(4,1,2,3)]
  return(nms_tab)
})

top5

saveRDS(top5, "./data/analysis/example_disease_top5_tissues.RDS")

latex_table <- xtable(top5[[1]])
print(latex_table, include.rownames = FALSE)

latex_table <- xtable(top10[[2]])
print(latex_table, include.rownames = FALSE)

latex_table <- xtable(top10[[3]])
print(latex_table, include.rownames = FALSE)

latex_table <- xtable(top10[[5]])
print(latex_table, include.rownames = FALSE)

latex_table <- xtable(top10[[6]])
print(latex_table, include.rownames = FALSE)
