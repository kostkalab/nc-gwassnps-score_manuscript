library(GenomicRanges)
library(seqminer)
library(rtracklayer)
library(ggplot2)

phegen.file  = "./external_data/phegen_dataset/PheGenI_Association_full.tab"
FILE_CHAIN   ="./external_data/ucsc/hg38ToHg19.over.chain"
# IPD          = "./external_data/1000_genomes/phase_3/"

# intermediate data read in
tss_19 = readRDS("./data/tss_19.RDS")
tss_38 = readRDS("./data/tss_38.RDS")

tss_382 = readRDS("~/project/variant-scores-2/data/genetic/tss_38.RDS")

phegen = read.table(phegen.file, 
                    sep = "\t", quote = "", comment.char = "", header = T, stringsAsFactors = F)


seqnames           = paste("chr",phegen$Chromosome,sep="")
seqnames[seqnames == "chr23"] = "chrX"
seqnames[seqnames == "chr24"] = "chrY"
seqnames[seqnames == "chr26"] = "chrM"
starts             = as.numeric(phegen$Location) + 1
phegen.38          = GRanges(seqnames,IRanges(start = starts, width = 1))
mcols(phegen.38 )  = phegen
coding             = c("missense", "cds-synon", "frameshift", "cds-indel", "STOP-GAIN") # fiter for noncoding
phegen.38          = phegen.38[!(phegen.38$Context %in% coding)]



# get nearest distance for ARB data
phegen.dist           = distanceToNearest(phegen.38, tss_38) # ARB data GR38
phegen.38$dist        = mcols(phegen.dist)[,1]

# map back to hg19
chain                = import.chain(FILE_CHAIN)
phegen.hg19          = unlist(GRangesList(liftOver(phegen.38,chain)))
names(phegen.hg19)   = paste0("rs", phegen.hg19$SNP.rs)

saveRDS(phegen.hg19, "./data/phegen_processed_hg19.RDS")
