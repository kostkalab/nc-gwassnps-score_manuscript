library(GenomicRanges)
library(dplyr)
library(tibble)

# read in phegen and GWAS and mapping information
catalog_all                 = readRDS("./data/GWAS_catalog_noncoding_hg19_uniquetrait.RDS")
phegen_all                  = readRDS("./data/phegen_processed_hg19.RDS")
mapped_trait_phegen_catalog = readRDS("./data/mesh_efo_mapping_45.RDS")
ARB_processed_hg19          = readRDS("./data/ARB_processed_hg19.RDS")


# make a combined object for all phegen and the catalog SNPs. for each trait, only keep unique snps
phenotypes = unique(mapped_trait_phegen_catalog$label)# only 41 traits

# as Dennis's instruction
# note on Jan 27: To get test and train snps (not considering neighbor snps), i can either get it overall or by disease. 
# i tried both ways, and it gives me back the same result. but since one snp can be annotated to more than one disease, it is easier to get them by disease. 
# Moreover, to move nearby snps between training and test set, i combine them together.- to make sure the test set is brand new: not only not seen by 
# DIVAN in the same disease, but overall not seen by DIVAN.


catalog_41 = lapply(phenotypes, function(x){
  print(x[1])
  efoname     = mapped_trait_phegen_catalog[mapped_trait_phegen_catalog$label %in% x, "mapped_label"] # get the efoname
  catalog     = catalog_all[catalog_all$UNIQUE.TRAIT %in% efoname] # extract GWAS catalog data
  
  catalog_new = catalog[!catalog %in% ARB_processed_hg19] # not in ARB
  catalog_new = catalog_new[catalog_new$DATE.ADDED.TO.CATALOG > "2016-05-28"] # after 2016-05-28
  
  names(catalog) = catalog$SNPS
  phegen  = phegen_all[phegen_all$Trait %in% x] # extract phegen data
  gr = c(granges(catalog), granges(phegen)) # combine phegen and the catalog
  gr = gr[!duplicated(gr)]# only keep unique SNPs
  if (length(gr) > 0) {    
    gr$trait = x
    # get a column indicating whether it's new or old
    gr$divan_new = 0
    if (length(catalog_new) > 0) {  gr[gr %in% catalog_new]$divan_new = 1}
  }
  return(gr)
})


ori = do.call("c",catalog_41)
a1 = ori[ori$divan_new == 1]
a0 = ori[ori$divan_new == 0]
dist = distanceToNearest(a1, a0)
dist_d = mcols(dist)[,"distance"]
a2 = a1[dist_d > 1000]
a0 = c(a0, a1[!dist_d > 1000]) # give back to a0
dist = distanceToNearest(a0, a2)
dist_d = mcols(dist)[,"distance"]
a0 = a0[dist_d > 1000] # remove ones in a0
a0$divan_new = 0

catalog_removed = c(a0,a2)

# look at disease-associated snps if they are next to each other. 
# remove snps in test snps (new snps) if there the distance is less than 1000
# catalog_41_removed = sapply(catalog_41, function(x){
#   a1 = x[x$divan_new == 1]
#   a0 = x[x$divan_new == 0]
#   if (length(a1) == 0| length(a0) == 0){
#     return(x)
#   }else{
#     dist = distanceToNearest(a1, a0)
#     dist_d = mcols(dist)[,"distance"]
#     a2 = a1[dist_d > 1000]
#     a0 = c(a0, a1[!dist_d > 1000])
#     dist = distanceToNearest(a0, a2)
#     dist_d = mcols(dist)[,"distance"]
#     a0 = a0[dist_d > 1000]
#     a0$divan_new = 0
#     print(c(length(a1)-length(a2), length(a1)))
#     return(c(a0, a2))
#   }
# })
# 
saveRDS(catalog_removed, "./data/divan_41.RDS")


# write out
cata_unlist = catalog_41 # i did this before i removed the nearby snps. keep this for now since having more snps doesn't hurt. 
names(cata_unlist) = NULL
cata_unlist = do.call("c", cata_unlist)
chr = gsub("chr","",seqnames(cata_unlist))
chr = gsub("X", "23", chr)

g.coord = paste0(chr,":",start(cata_unlist))
g.coord = c("lead_snp",unique(g.coord))

write.table(g.coord, 
            file = "./snpsnap/input_divan_41.txt", 
            sep = "\n", row.names = F, col.names = F, quote = F)


