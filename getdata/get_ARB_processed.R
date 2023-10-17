library(GenomicFeatures)
library(GenomicRanges)
library(pegas)
library(biomaRt)
library(VariantAnnotation)


# data FILE
ARB_FILE     = "./external_data/ARB_dataset/ARB_GWAS_SNPS_20190115.TAB"
FILE_CHAIN   ="./external_data/ucsc/hg38ToHg19.over.chain"
# IPD          = "/data/projects/annotation/1000_genomes/phase_3/"

# intermediate data read in
tss_19 = readRDS("./data/tss_19.RDS")
tss_38 = readRDS("./data/tss_38.RDS")




#------------ 
# GET ARB GWAS DATA
# read in ARB GWAS data
AW                 = read.table(ARB_FILE, sep = "\t", header = T, comment.char = "", quote = "", stringsAsFactors = F)
seqnames           = paste("chr",AW$Chr,sep="")
starts             = as.numeric(AW$Position)
ARB.SNPS.38        = GRanges(seqnames,IRanges(start = starts, width = 1))
mcols(ARB.SNPS.38) = AW
coding             = c("Missense", "Cds-synon", "Frameshift") # fiter for noncoding
ARB.SNPS.38        = ARB.SNPS.38[!(ARB.SNPS.38$Context %in% coding)]


# get nearest distance for ARB data
ARB.dist             = distanceToNearest(ARB.SNPS.38, tss_38) # ARB data GR38
ARB.SNPS.38$dist     = mcols(ARB.dist)[,1]

chain                = import.chain(FILE_CHAIN)
ARB.SNPS.hg19        = unlist(GRangesList(liftOver(ARB.SNPS.38,chain)))
names(ARB.SNPS.hg19) = ARB.SNPS.hg19$SNP 

## there are duplicated snps


# ARB.SNPS.hg19.dis    = disjoin(ARB.SNPS.hg19) # some SNPs are discovered by several studies, make them into 1
# ARB.SNPS.hg19        = ARB.SNPS.hg19[match(ARB.SNPS.hg19.dis, ARB.SNPS.hg19)] # 20269 SNPs left

#--------------------
# bonus
## get reference and alternative allele for ARB.SNPS.hg19
#  IPF = Sys.glob(paste(IPD,"*genotypes.vcf.gz",sep=""))
#  ARB_alt = lapply(IPF, function(x){
#    a = ARB.SNPS.hg19
#    seqlevelsStyle(a) = "NCBI"
#    chr = gsub("ALL.chr|.phase3", "", strsplit(basename(x), "_")[[1]][1])
#    a = a[seqnames(a) ==  chr]
#    svp <- ScanVcfParam(info=c("AF", "VT"), geno=NA, which=a) # only read in AF and VT
#    tab <- TabixFile(x)
#    vcf1 = readVcf(tab, "hg19", svp)
#    cat("readin",x, "\n")
#    vcf2 = vcf1[isSNV(vcf1)]
#    r = rowRanges(vcf2)
#    return(r)
#  })
#  ARB_alt = do.call("c", ARB_alt)
#  seqlevelsStyle(ARB_alt) = "UCSC"
#  ARB_alt = ARB_alt[!width(ARB_alt) > 1]
#  f = findOverlaps(ARB.SNPS.hg19, ARB_alt) # same number
#  ARB_alt[subjectHits(f)]
#  
#  ARB.SNPS.hg19$REF = NA_character_
#  ARB.SNPS.hg19$ALT = NA_character_
#  ARB.SNPS.hg19[match(ARB_alt[subjectHits(f)], ARB.SNPS.hg19)]$REF = ARB_alt[subjectHits(f)]$REF
#  ARB.SNPS.hg19[match(ARB_alt[subjectHits(f)], ARB.SNPS.hg19)]$ALT = ARB_alt[subjectHits(f)]$ALT


saveRDS(ARB.SNPS.hg19, file = "./data/ARB_processed_hg19.RDS")
