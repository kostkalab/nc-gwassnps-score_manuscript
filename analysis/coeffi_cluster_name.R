library(hasseDiagram)
library(ontologyIndex) # remember to use version 2.2 to have the get_OWL function
library(readr)
file = "./external_data/efo_ontology/3.24.0/efo_3_24_0.obo"
efo = get_ontology(file, propagate_relationships = "is_a", extract_tags = "minimal") # i'm pretty sure rols package use is_a not part_of
a1 = read_csv("./sup_data/sup_data_disease-terms.csv.gz")



coeffi_clusters <- readRDS("~/project/variant-scores_manuscript/data/analysis/coeffi_clusters.RDS")
a21 = coeffi_clusters[,c("disease", "clst", "clst_nm")] # choose the useful columns
colnames(a21) = c("phenotype", "cluster_id", "cluster_name") # rename columns

get_lower = function(x){
  y = tolower(gsub(" ", "_", x))
  return(y)
}

a1$term_name = get_lower(a1$term_name)
a21$efo = a1$term_efo_id[match(a21$phenotype, a1$term_name)]



get_matrix_with_names = function(ontology, terms){
  mtx = ontologyIndex::get_term_descendancy_matrix(ontology, terms)
  rownames(mtx) = efo$name[rownames(mtx)]
  colnames(mtx) = efo$name[colnames(mtx)]
  return(mtx)
}


exclude = c("disposition", "material property", "experimental factor", "disease" , 
            "disease by anatomical system", "disorder by anatomical region")
exclude = c("MONDO:0021199", "EFO:0000408", "EFO:0000001", "BFO:0000020", "BFO:0000016", 
            "MONDO:0024505")
get_frequency = function(cluster){
  print(cluster)
  acl = a21[a21$cluster_name == cluster,]
  matrix1 = get_matrix_with_names(efo, acl$efo) # get the matrix of who include who
  acl$efo
  
  fre = get_term_frequencies(efo, acl$efo, patch_missing = FALSE)
  fre = fre[!names(fre) %in% exclude]
  
  top_fre = sort(fre, decreasing = T)[1:10]
  df = data.frame(name = efo$name[names(top_fre)], term_id = names(top_fre),frac = top_fre)
  df$name = factor(df$name, levels = df$name)
  df$cluster = cluster
  return(df)
}

dfs = lapply(unique(a21$cluster_name), get_frequency)

dfs = do.call("rbind", dfs)
rownames(dfs) = NULL

colnames(dfs) = c("term", "term_id","term_frequency", "cluster_name")

dfs %>% dplyr::mutate(dplyr::across(term_frequency, format, digits = 4))  %>% 
  readr::write_csv(file="./sup_data/sup_data_cluster_term-frequency.csv.gz") 


# 
# nervous = a1[a1$term_efo_id %in% get_descendants(efo, "EFO:0000618"),]
# mtx = get_matrix_with_names(efo,nervous$term_efo_id)
# hasse(mtx)
