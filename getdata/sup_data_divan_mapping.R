library(readr)
# DIVAN mapping
mapped_trait_phegen_catalog = readRDS("./data/mesh_efo_mapping_45.RDS")
colnames(mapped_trait_phegen_catalog) = c("MESH_ID", "MESH_LABEL", "EFO_ID", "EFO_LABEL", "distance")

readr::write_csv(mapped_trait_phegen_catalog[1:4], file="./sup_data/sup_data_mapping-efo-mesh.csv.gz")

