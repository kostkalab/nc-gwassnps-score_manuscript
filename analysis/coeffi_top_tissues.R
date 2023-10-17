library(reshape2)

# get coefficients and epigenome annotation
#------------------------------------------

coeff              <- readRDS("./data/analysis/coefficients.RDS")
#- get mean of coefficients
top_mean           <- sapply(coeff, function(x) {x[,1]})
rownames(top_mean) <- rownames(coeff[[1]])
coeffi             <-  top_mean[-c(1,2),] # remove rowmeans and intercept, only concentrate on tissues

# retrieve ANATOMY, TYPE, and GROUP for 127 epigenomes. 
library(erma) # https://bioconductor.org/packages/devel/bioc/vignettes/erma/inst/doc/erma.html
meta          <- mapmeta()
met           <- as.data.frame(meta[,c("Epigenome.ID..EID.", "GROUP","Standardized.Epigenome.name", "ANATOMY", "TYPE")])
rownames(met) <- met[,1]
met           <- met[rownames(coeffi),]

tt = readRDS("./data/analysis/coeffi_clusters.RDS")

######### start to get top 5 tissues
coeffi_norm = apply(coeffi,2,function(x) {x=x-min(x); x=  x/quantile(x,.95); x[x>1]=1; return(x)}) # QQ: kind of normalize it; anything larger than 1 is 1
co_scale = coeffi_norm


# get the p value for each tissue for each cluster
gpv <- function(clind) {
  rs = apply(co_scale,1, function(rw) {
    rw2 = as.data.frame(cbind(rw, group = tt$clst)) # add group name
    gm = tapply(rw2$rw, rw2$group, mean)# get the group mean
    best = names(gm[names(gm) != clind][which.max(gm[names(gm) != clind])]) # exclude the cluster
    wilcox.test(x = rw[tt$clst == clind], y = rw[tt$clst == best], alternative = "greater")
  })
  pv = lapply(rs,function(x) x$p.value)
  return(unlist(pv))
} # get the p-value for each cluster

# gpv <- function(clind) {
#   rs = apply(coeffi,1, function(rw) {
#     wilcox.test(x = rw[tt$clst == clind], y = rw[tt$clst != clind], alternative = "greater")
#   })
#   pv = lapply(rs,function(x) x$p.value)
#   return(unlist(pv))
# } # get the p-value for each cluster


wpvs = sapply(1:7,gpv) # get the pvalue for each cluster


# find top 5 tissues
top_5 = apply(wpvs, 2, function(x){
  ord = order(x)[1:5]
  tiss_5 = rownames(wpvs[ord,])
  return(tiss_5)
})
df_top5 = as.data.frame(t(top_5)) # make it a data.frame
df_top5$cluster = rownames(df_top5)
clusterdf = melt(df_top5, id.vars = "cluster")
names(clusterdf) = c("cluster_id", "rank", "tissue")
clusterdf = clusterdf[order(clusterdf$cluster),]

# expand it. 


## df_cluster_nm = data.frame(clustid = 1:7,
##                            longname = c("immune-1", "others", "cardiovascular, others", "immune-2", 
##                                         "mental or behavioural disorder", "digestive system disease, cancer", "skin cancer"), 
##                            shortname = c("immune-1", "others", "cardio,others", "immune-2", 
##                                          "mental", "digest,cancer", "skin cancer"))

x1 = sapply(1:7, function(x){
  y = tt[tt$clst == x,][1,]
  y1 = y[,c("clst", "clst_nm", "clst_shortnm")]
  colnames(y1) = c("clustid","longname","shortname")
  return(y1)
})
df_cluster_nm = as.data.frame(t(x1))

clusterdf$cluster_nm = df_cluster_nm[match(clusterdf$cluster_id, df_cluster_nm$clustid),"longname"]

## prepare to write out
m = match(clusterdf$tissue, met$Epigenome.ID..EID.)
clusterdf$tissue_name = met[m,"Standardized.Epigenome.name"]
clusterdf$tissue_anatomy = met[m, "ANATOMY"]

# rearrange the order of columns
clusterdf = clusterdf[,c("cluster_id", "cluster_nm", "rank", "tissue", "tissue_name", "tissue_anatomy")]

saveRDS(clusterdf, file = "./data/analysis/coeffi_top5tissues.RDS")

