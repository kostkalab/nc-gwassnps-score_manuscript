library(GenomicRanges)
library(readr)

# read in objects
snpsnap_catalog_ct100_aggre <- readRDS("./data/catalog_cutoff100_snpsnap.RDS")
snpsnap_catalog_ct100_aggre = snpsnap_catalog_ct100_aggre[order(snpsnap_catalog_ct100_aggre$phenotype)] # order by phenotype
snpsnap_tss_catalog_ct100_aggre <- readRDS("./data/catalog_cutoff100_snpsnap_tss.RDS")
snpsnap_tss_catalog_ct100_aggre = snpsnap_tss_catalog_ct100_aggre[order(snpsnap_tss_catalog_ct100_aggre$phenotype)]
tss_catalog_ct100_aggre <- readRDS("./data/catalog_cutoff100_tss.RDS")
tss_catalog_ct100_aggre = tss_catalog_ct100_aggre[order(tss_catalog_ct100_aggre$phenotype)]
random_catalog_ct100_aggre <- readRDS("./data/catalog_cutoff100_random.RDS")
random_catalog_ct100_aggre = random_catalog_ct100_aggre[order(random_catalog_ct100_aggre$phenotype)]

ld_marked = readRDS("./data/catalog_cutoff100_snpsnap_ld_marked.RDS")
ld_marked = ld_marked[order(ld_marked$phenotype)]

# generate suppl2: disease associated SNVs - LD block marked 
disease             = ld_marked[ld_marked$type == 1]
df_disease = data.frame(SNV_ID = disease$snpID, 
                        rsID=disease$rsID,
                        phenotype = disease$phenotype, 
                        hg19_chromosome = seqnames(disease), 
                        hg19_location = start(disease), 
                        LD_block_cluster = disease$cluster, 
                        representive_SNP = disease$ld_block_picked, 
                        stringsAsFactors = F)




# generate suppl3: control SNVs
control_snpsnap     = snpsnap_catalog_ct100_aggre[snpsnap_catalog_ct100_aggre$type == 0]
control_snpsnap_tss = snpsnap_tss_catalog_ct100_aggre[snpsnap_tss_catalog_ct100_aggre$type == 0]
control_tss         = tss_catalog_ct100_aggre[tss_catalog_ct100_aggre$type == 0]
control_random      = random_catalog_ct100_aggre[random_catalog_ct100_aggre$type == 0]

control_snpsnap$matched_method     = "snpsnap"
control_snpsnap_tss$matched_method = "snpsnap_tss"
control_tss$matched_method         = "tss"
control_random$matched_method      = "random"          

control = c(control_snpsnap, control_snpsnap_tss, control_tss, control_random)

df_control = data.frame(SNV_ID = control$snpID, 
                        rsID = control$rsID,
                        phenotype = control$phenotype, 
                        hg19_chromosome = seqnames(control), 
                        hg19_location = start(control), 
                        matched_method = control$matched_method,
                        disease_SNV_ID = control$input_snp, 
                        stringsAsFactors = F)


readr::write_csv(df_disease, "./sup_data/sup_data_disease-snvs.csv.gz")
readr::write_csv(df_control, "./sup_data/sup_data_control-snvs.csv.gz")
