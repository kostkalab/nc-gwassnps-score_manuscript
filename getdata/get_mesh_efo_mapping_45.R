phegen_45 <- readRDS("~/project/variant-scores-2/data/phegen_project/genetic/phegen_45.RDS")
phegen_45$phenotype = droplevels(phegen_45$phenotype)

library(GenomicRanges)

ESM = read.csv("./external_data/DIVAN_paper/13059_2016_1112_MOESM1_ESM-3_TableS2.csv")
mesh.names = unique(droplevels(phegen_45$phenotype))

phegen_processed_hg19 <- readRDS("~/projects/variant-scores-2/data/phegen_project/genetic/phegen_processed_hg19.RDS")

map_library = read.table("./external_data/oxo_mapping/oxo_mesh_efo_mapping_distance3.txt",
                         sep = ",", header = T)


mapped_EFO = map_library[map_library$label %in% mesh.names, c(1,2,3,4,7)]
mapped_EFO1 = mapped_EFO[mapped_EFO$distance == 1,] # some mesh terms are mapped with 2 efo terms
# mapped_EFO1$label = droplevels(mapped_EFO1$label)
# mapped_EFO1$mapped_label = droplevels(mapped_EFO1$mapped_label)

saveRDS(mapped_EFO1, "./data/mesh_efo_mapping_45.RDS")
