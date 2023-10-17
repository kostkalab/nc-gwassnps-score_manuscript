library(readr)

tt = readRDS("~/project/variant-scores-2/analysis/GWAS_catalog_3.0/aggre/result_data/111_diseases_cluster.RDS")
tt = readRDS("./data/analysis/coeffi_clusters.RDS")
suppl21 = tt[,c("disease", "clst", "clst_nm")]
colnames(suppl21) = c("phenotype", "cluster_id", "cluster_name")
readr::write_csv(suppl21, file="./sup_data/sup_data_cluster-id-name.csv.gz") 


clusterdf = readRDS("./data/analysis/coeffi_top5tissues.RDS")

readr::write_csv(clusterdf[,c(1:2,4:6)], file="./sup_data/sup_data_top-five-tissues.csv.gz") 



