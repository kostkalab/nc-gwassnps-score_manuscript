library(hasseDiagram)
library(ontologyIndex) # remember to use version 2.2 to have the get_OWL function
library(readr)
file = "./external_data/efo_ontology/3.24.0/efo_3_24_0.obo"
efo = get_ontology(file, propagate_relationships = "is_a", extract_tags = "minimal") # i'm pretty sure rols package use is_a not part_of


coeffi_clusters <- readRDS("./data/analysis/coeffi_clusters.RDS")

a1 = read_csv("./sup_data/sup_data_disease-terms.csv.gz")

a21 = coeffi_clusters[,c("disease", "clst", "clst_nm")] # choose the useful columns
colnames(a21) = c("phenotype", "cluster_id", "cluster_name") # rename columns

get_lower = function(x){
  y = tolower(gsub(" ", "_", x))
  return(y)
}


a1$term_name = get_lower(a1$term_name)
a21$efo = a1$term_efo_id[match(a21$phenotype, a1$term_name)]

get_matrix_with_names = function(ontology, terms){
  mtx = get_term_descendancy_matrix(ontology, terms)
  rownames(mtx) = get_lower(efo$name[rownames(mtx)])
  colnames(mtx) = get_lower(efo$name[colnames(mtx)])
  return(mtx)
}

pdf(file="./plot/cluster_hesse_new.pdf", width = 15, height = 7)
plot_hasse = function(cluster_name){
  a21_cluster = a21[a21$cluster_name == cluster_name,]
  add_terms = add_terms_list[[cluster_name]]
  efoid = efo$id[efo$name %in% add_terms]
  all_id = unique(c(a21_cluster$efo, efoid))
  mtx = get_matrix_with_names(efo, all_id)
  label = rownames(mtx)
  label[!label %in% a21$phenotype] = paste(label[!label %in% a21$phenotype],"**", sep = "")
  label[label %in% a21$phenotype & !label %in% a21_cluster$phenotype] = 
    paste(label[label %in% a21$phenotype & !label %in% a21_cluster$phenotype], "*",sep = "")
  p = hasse(mtx,label = label, list(cluster = T))
}


cluster_names = unique(a21$cluster_name)
add_terms_list = list(`immune system disease` = c("immune system disease", "asthma"), 
                      heterogenous = c("nervous system disease", "skeletal system disease",
                                 "lung disease",
                                 "metabolic disease","abdominal and pelvic region disorder"),
                      `cardiovascular disease/others` = c("cardiovascular disease", "nervous system disease"), 
                      `immune system disease/autoimmune disease` = NULL, 
                      `mental or behavioural disorder` = "nervous system disease", 
                      `digestive system disease/cancer` = c("cancer", "digestive system disease"),
                      `skin cancer` = c("skin cancer"))

lapply(cluster_names, plot_hasse)
dev.off()
