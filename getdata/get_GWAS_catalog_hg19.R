library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)

FILE_CATALOG="~/project/variant-scores_manuscript/external_data/gwas_catalog/gwas_catalog_v1.0.2-associations_e100_r2020-12-02.tsv.gz"
FILE_CHAIN="~/project/variant-scores_manuscript/external_data/ucsc/hg38ToHg19.over.chain"

GW = read.csv(FILE_CATALOG,sep="\t",stringsAsFactors=FALSE, quote = "")

# filter 1: only keep SNVs
#- make a GRanges object from it / remove multi-location entries (cheap)
seqnames    = paste("chr",GW$CHR_ID,sep="")
starts      = as.numeric(GW$CHR_POS)
#- remove multi-loctation entries
id          = !is.na(starts)
seqnames    = seqnames[id]  # delete id which is NA
starts      = starts[id]
ends        = starts
SNPS        = GRanges(seqnames,IRanges(starts,ends))
mcols(SNPS) = GW[id,]

#- move SNPS from hg38 to hg19
chain     = import.chain(FILE_CHAIN)
SNPS.hg19 = unlist(GRangesList(liftOver(SNPS,chain))) # all single location entries of GWAS catalog with hg19

saveRDS(SNPS.hg19, "~/project/variant-scores_manuscript/data/GWAS_catalog_hg19.RDS")

# filter 2: only keep non-coding SNVs
SNPS.hg19 = readRDS("~/project/variant-scores_manuscript/data/GWAS_catalog_hg19.RDS")
exon_term = c("missense_variant", "synonymous_variant", "frameshift_variant", 
              "stop_gained", "stop_lost", "start_lost", 
              "inframe_deletion", "inframe_insertion",
              "coding_sequence_variant", "protein_altering_variant") # double check this


# keep variants that locate on noncoding regions
noncoding = readRDS("~/project/variant-scores_manuscript/data/noncoding_hg19.RDS")
f = findOverlaps(SNPS.hg19, noncoding)
SNPS_noncoding = SNPS.hg19[queryHits(f)]

# keep the variants with context not in exon terms. 
SNPS.hg19_noncoding =  SNPS_noncoding[!SNPS_noncoding$CONTEXT %in% exon_term]
length(SNPS.hg19_noncoding) / length(SNPS.hg19)

saveRDS(SNPS.hg19_noncoding, "~/project/variant-scores_manuscript/data/GWAS_catalog_noncoding_hg19.RDS")
