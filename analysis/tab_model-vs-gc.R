a1 = read_csv("./sup_data/sup_data_beta-model-similarity-dhs.csv.gz")
a2 = read_csv("./sup_data/sup_data_genetic-correlation.csv.gz")

a3 = merge(a2,a1, by=c("disease_1","disease_2")) # NA's match

# set the threshold to be 0.5
## a4 = a3 %>% dplyr::filter(genetic_correlation > 0.5 | model_similarity > 0.5) 
## 
## a4$group = NA
## a4 %>% 
##   mutate(group = replace(group, genetic_correlation > 0.25 & model_similarity > 0.25, "B")) %>%
##   mutate(group = replace(group, genetic_correlation < 0.25 & model_similarity > 0.25, "C")) %>%
##   mutate(group = replace(group, genetic_correlation > 0.25 & model_similarity < 0.25, "D")) %>%
##   arrange(group)

quantile(a3$genetic_correlation, c(0.9)) # gc pick 0.20
quantile(a3$model_similarity, c(0.9)) # model_sim pick 0.25

a3$group = NA
a3 = a3 %>% 
  mutate(group = replace(group, genetic_correlation > 0.20 & model_similarity > 0.25, "B")) %>%
  mutate(group = replace(group, genetic_correlation < 0.20 & model_similarity > 0.25, "C")) %>%
  mutate(group = replace(group, genetic_correlation > 0.20 & model_similarity < 0.25, "D")) %>%
  arrange(group)
  
a3[!is.na(a3$group),]
  
b1 = a3 %>% dplyr::filter(genetic_correlation > 0.20 & model_similarity > 0.25)
b1$diff = b1$genetic_correlation + b1$model_similarity
b1 = b1[order(b1$diff, decreasing = T),]
b3 = b1[1:8,]
c1 = a3 %>% dplyr::filter(genetic_correlation > 0.20 & model_similarity < 0.25)
c1$diff = c1$genetic_correlation-c1$model_similarity
c1 = c1[order(c1$diff, decreasing = T),]
c3 = c1[1:8,]
d1 = a3 %>% dplyr::filter(genetic_correlation < 0.20 & model_similarity > 0.25)
d1$diff = d1$model_similarity-d1$genetic_correlation
d1 = d1[order(d1$diff, decreasing = T),]
d3 = d1[1:8,]

a4 = rbind(b3,c3,d3)

tab_examp = apply(a4, 1, function(x) x[1:2])
id2 = apply(tab_examp, 2,find_exmp)

a4[,1:5] %>% dplyr::mutate(dplyr::across(c(3), format, digits = 2))  %>% 
  dplyr::mutate(dplyr::across(c(4), format, digits = 2))  %>% 
  readr::write_csv(file = "./sup_data/tab_model-vs-gc.csv") 

firstup = function(x){
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


a4$disease_1 = gsub("_", " ", firstup(a4$disease_1))
a4$disease_2 = gsub("_", " ", firstup(a4$disease_2))
rownames(a4) = NULL
print(xtable(a4[,1:5], type = "latex"), include.rownames=FALSE)

