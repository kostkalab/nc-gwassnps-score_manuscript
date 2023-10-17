library(dplyr)
library(tidyr)

#- notes by Dennis
### notes by QQ

gwas_atls <- readr::read_tsv("./external_data/gwas_atlas/gwasATLAS_v20191115.txt.gz", guess_max=10000)


gwas_atls %>% head(n=3)

gwas_atls_gc <- readr::read_tsv("./external_data/gwas_atlas/gwasATLAS_v20191115_GC.txt.gz", guess_max=10000)

tmp <- gwas_atls %>% semi_join(gwas_atls_gc, by = c("id"="id1"))
tmp <- rbind(gwas_atls %>% filter(id==1), tmp) #- will have missed the first one 
gwas_atls_fil <- tmp ### not all the gwas atals id have gc, find the ones with gc

qq_pt  <- readr::read_csv("./sup_data/sup_data_disease-terms.csv.gz")
qq_snp <- readr::read_csv("./sup_data/sup_data_disease-snvs.csv.gz",
                          col_types = readr::cols( SNV_ID = readr::col_character(),
                                                   rsID  = readr::col_character(),
                                                   phenotype = readr::col_character(),
                                                   hg19_chromosome = readr::col_character(),
                                                   hg19_location = readr::col_integer(),
                                                   LD_block_cluster = readr::col_integer(),
                                                   representive_SNP = readr::col_integer()))

qq_all <- full_join(qq_snp,qq_pt, by = c("phenotype" = "term_name"))

qq_all %>% head(n=3)

#- to lowercase and strip punctuation, collapse blank(s)
strng_tls <- function(strng){
  #============================
  gsub("[[:punct:][:blank:]]+", " ", tolower(strng))
}

map_efo <- function(efo_id, gwas_atls){
  #=======================================
  
  #- get synonyms
  trm = rols::term("efo",efo_id)
  syn = rols::termSynonym(trm)
  syn = strng_tls(syn)
  syn = c(strng_tls(rols::termLabel(trm)),syn)
  
  #- perform matching
  ### agrep https://docs.tibco.com/pub/enterprise-runtime-for-R/4.4.0/doc/html/Language_Reference/base/agrep.html
  ind1 = unlist(sapply(paste("^",syn,"$",sep=""), agrep,
                       strng_tls(gwas_atls$uniqTrait), max.distance=list(substitutions=1), fixed=FALSE)) ### allow only one substution
  ind2 = unlist(sapply(paste("^",syn,"$",sep=""), agrep,
                       strng_tls(gwas_atls$Trait), max.distance=list(substitutions=1), fixed=FALSE))
  inds = sort(unique(c(ind1,ind2)))
  
  if(length(inds) ==0) {
    return(data.frame(efo=efo_id, atl_ids = NA))        
  } #- otherwise proper result 
  return(data.frame(efo = efo_id, atl_ids = gwas_atls$id[inds]))
  
}

#- map and make tibble
mtch <- sapply(sort(unique(qq_pt$term_efo_id)), map_efo, gwas_atls_fil) 
mtch <- as_tibble(plyr::adply(mtch,2,function(x) data.frame(x))) %>% select(efo_term_id = efo,gwas_atls_id = atl_ids)
### table(is.na(mtch$gwas_atls_id))
### among 111 diseases, we found 35 diseases in GWAS atlas. (the same number as my code, checked!)
### one efo term can be mapped to more than one atlas terms

#- add readable names                              
mtcha <- left_join(mtch, qq_pt %>% select(term_efo_id,term_name), by=c("efo_term_id" = "term_efo_id"))
mtcha <- left_join(mtcha, gwas_atls_fil %>% select(id,uniqTrait), by=c("gwas_atls_id" = "id"))
mtcha <- mtcha %>% mutate(term_name = strng_tls(term_name), uniqTrait = strng_tls(uniqTrait))                     

mtcha %>% head(n=22)

### only keep the mapped ones: mtchb
mtchb <- left_join(mtcha, gwas_atls_fil %>% select(id,gwas_atls_N=N), by=c("gwas_atls_id" = "id"))
mtchb <- mtchb %>% filter(!is.na(gwas_atls_id))
mtchb %>% head 

### add the weights of each study: mtchc
tmp   <- mtchb %>% group_by(efo_term_id) %>% summarize(nids = n(), gwas_atls_Ntotal = sum(gwas_atls_N))
mtchc <- left_join(mtchb,tmp,by="efo_term_id") %>% mutate(wght = gwas_atls_N/gwas_atls_Ntotal) %>%
  select(-c(nids, gwas_atls_Ntotal))
mtchc %>% head()

### check if the weight is 1
tmp <- mtchc %>% group_by(efo_term_id) %>% summarize(sw = sum(wght)) 
tmp %>% head()
print(all(tmp$sw -1 < sqrt(.Machine$longdouble.eps))) ### may not be important

tmp <- t(combn(sort(mtchc$efo_term_id %>% unique()),2)) %>% as.data.frame() %>% tibble() ### combination of two efo terms
colnames(tmp) = c("efo_term_id_1","efo_term_id_2")
tmp$efo_term_id_1 = as.character(tmp$efo_term_id_1) ### QQ add: factor doesn't work, so i changed it to charactor
tmp$efo_term_id_2 = as.character(tmp$efo_term_id_2) ### same above
tmp <- tmp %>% filter(efo_term_id_1 != efo_term_id_2) ### stays the same number actually...
tmp %>% head()


tmp2 <- left_join(tmp,  mtchc %>% select(efo_term_id_1 = efo_term_id, wght_1 = wght, gwas_atls_id_1= gwas_atls_id),  by=c("efo_term_id_1" = "efo_term_id_1"))
tmp3 <- left_join(tmp2, mtchc %>% select(efo_term_id_2 = efo_term_id, wght_2 = wght, gwas_atls_id_2 = gwas_atls_id), by=c("efo_term_id_2" = "efo_term_id_2"))
tmp4 <- tmp3 %>% rowwise() %>% mutate(gwas_atls_id_sml = min(gwas_atls_id_1,gwas_atls_id_2), gwas_atls_id_mx = max(gwas_atls_id_1,gwas_atls_id_2))
### note for tmp4: in gwasatlas gc, id1 is always smaller than id2. so make two columns, with the first one's id always smaller than the second one
tmp5 <- left_join(tmp4, gwas_atls_gc %>% select(id1, id2, rg, se, p), by=c("gwas_atls_id_sml" = "id2", "gwas_atls_id_mx" = "id1"))

efo_gc_big <- tmp5; #rm(tmp,tmp2,tmp3,tmp4,tmp5)
efo_gc_big %>% head(n=5)

tmp <- efo_gc_big %>% group_by(efo_term_id_1, efo_term_id_2) %>% summarize(sw = sum(wght_1*wght_2))
tmp %>% head(n=5) ### make sure the weights add up to one
print(all(tmp$sw-1.0 < sqrt(.Machine$longdouble.eps))) #- comparing to ~ 3.3E-10

#- note that NA's from gwas_atls_gc propagated into efo_gc_big...

sum_and_renorm <- function(wts,rg){
  #==================================
  ind = !is.na(rg)
  wts = wts[ind]
  rg  = rg[ind]
  wts = wts/sum(wts)
  return(sum(wts*rg))
}

efo_gc <- efo_gc_big %>% group_by(efo_term_id_1,efo_term_id_2) %>% summarize(rg = sum_and_renorm(wght_1*wght_2,rg))
efo_gc <- left_join(efo_gc, qq_pt %>% select(term_efo_id, term_name_id_1 = term_name), by=c("efo_term_id_1" = "term_efo_id"))
efo_gc <- left_join(efo_gc, qq_pt %>% select(term_efo_id, term_name_id_2 = term_name), by=c("efo_term_id_2" = "term_efo_id"))
efo_gc %>% head()

efo_gc %>% arrange(-rg) %>% head(n=22) ### look at the top 22 pairs of genetic correlation
barplot(sort(efo_gc$rg))

saveRDS(efo_gc, "./data/analysis/efo_gc.RDS")
## readr::write_csv(efo_gc_big %>% select(-gwas_atls_id_sml, -gwas_atls_id_mx)  %>% 
##                    mutate(across(where(is.numeric), round, 6)), 
##                  file="./data/analysis/efo_gc_big.csv")
## readr::write_csv(efo_gc     %>% mutate(across(where(is.numeric), round, 6)), 
##                 file="./data/analysis/efo_gc.csv")


