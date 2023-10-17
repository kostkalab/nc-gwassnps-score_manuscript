library(readr)
library(dplyr)

smat = readRDS("./data/analysis/coeffi_weighted_correlations.RDS")
sdat           <- data.frame(t(combn(rownames(smat),2)), stringsAsFactors = F) ### don't know why it doesn't work for me. have to change to charactor to work
sdat$sim       <- smat[cbind(sdat[,1],sdat[,2])] ### Nice way to do it!!
sdat           <- dplyr::tibble(sdat)
colnames(sdat) <- c("term_1","term_2","sim")
head(sdat)

gdat <- readRDS("./data/analysis/efo_gc.RDS")



sdat <- sdat %>% rowwise() %>% mutate(id_1 = min(term_1,term_2), id_2 = max(term_1,term_2))
sdat2 = sdat %>% dplyr::select(id_1, id_2, sim) %>% 
  dplyr::rename(disease_1=id_1, disease_2=id_2, model_similarity=sim)

head(sdat2)
readr::write_csv(sdat2, file=("./sup_data/sup_data_beta-model-similarity-dhs.csv.gz")) 



strng_tls <- function(strng){
  #============================
  gsub("[[:punct:][:blank:]]+", " ", tolower(strng)) %>% stringr::str_replace_all("[[:space:]]+","_")
}

gdat <- gdat %>% dplyr::mutate(term_name_id_1 = term_name_id_1 %>% tolower %>%  stringr::str_replace_all("[[:space:]]+","_"))
gdat <- gdat %>% dplyr::mutate(term_name_id_2 = term_name_id_2 %>% tolower %>%  stringr::str_replace_all("[[:space:]]+","_"))

gdat %>% head()
gdat <- gdat %>% dplyr::rowwise() %>% dplyr::mutate(id_1 = min(term_name_id_1,term_name_id_2), id_2 = max(term_name_id_1,term_name_id_2))
gdat %>% head() 



gdat2 = gdat %>% ungroup() %>% dplyr::select(id_1, id_2, rg) %>% 
  dplyr::rename(disease_1=id_1, disease_2=id_2, genetic_correlation=rg)
gdat2 %>% head() 
readr::write_csv(gdat2, file="./sup_data/sup_data_genetic-correlation.csv.gz") 
