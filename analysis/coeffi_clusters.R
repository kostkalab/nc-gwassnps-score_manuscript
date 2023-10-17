library(uwot)
library(dplyr)

## in this script, we assume the cluster name is already known (which is not technically true). 
## Please see the other file, "./analysis/coeffi_cluster_name.R" to find out the detail on how we give the name to the clusters.

tmp = readRDS("./data/analysis/coeffi_weighted_correlations.RDS")

### do umap
set.seed(123432)
tu0 <- umap(tmp,ret_model=TRUE,pca_center = FALSE, bandwidth = 1.5) # do umap
rownames(tu0$embedding) <- colnames(tmp)

### cluster into 7 clusters
h = hclust(dist(tu0$embedding))
clus <- cutree(h, 7)  # seven clusters


tt <- as.data.frame(tu0$embedding)
colnames(tt) = c("d1","d2")
tt$disease = rownames(tt) # get disease names

tt$clst = as.factor(clus)
rownames(tt) = NULL # remove the rownames


df_cluster_nm = data.frame(clustid = 1:7,
                           longname = c("immune system disease", "heterogenous", "cardiovascular disease/others", 
                                        "immune system disease/autoimmune disease", 
                                        "mental or behavioural disorder", "digestive system disease/cancer", "skin cancer"), 
                           shortname = c("immune", "heterogenous", "cardio/others", "immune/autoimmune", 
                                         "mental", "digest/cancer", "skin cancer"))


# add a column to represent the name of the clusters
match(tt$clst, df_cluster_nm$clustid)
tt$clst_nm = df_cluster_nm$longname[tt$clst] # the match happens to be the same as the cluster id, so use it for simplicity
tt$clst_shortnm = df_cluster_nm$shortname[tt$clst]

clus2 = cutree(h, 2)  # two clusters
tt$clus2 = factor(clus2)
tt$clst2_nm = NA
tt = tt %>%
  mutate(clst2_nm=replace(clst2_nm, clus2==1, "immune")) %>%
  mutate(clst2_nm=replace(clst2_nm, clus2==2, "others"))

saveRDS(tt, "./data/analysis/coeffi_clusters.RDS")
