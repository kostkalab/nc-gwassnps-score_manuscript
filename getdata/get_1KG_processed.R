library(GenomicFeatures)
library(GenomicRanges)
library(pegas)
library(biomaRt)
library(VariantAnnotation)
library(GenomeInfoDb)



# data FILE
IPD          = "~/project/variant-scores_manuscript/external_data/1000_genomes/phase_3/"

# intermediate data read in
tss_19 = readRDS("~/project/variant-scores_manuscript/data/tss_19.RDS")
tss_38 = readRDS("~/project/variant-scores_manuscript/data/tss_38.RDS")
# ARB.SNPS.hg19 = readRDS("~/projects/variant-scores-2/data/genetic/ARB_processed_hg19.RDS")


#----------------------
# GET 1000 GENOME DATA
# read in vcf file, and filter for SNP and allele frequency
IPF = Sys.glob(paste(IPD,"*genotypes.vcf.gz",sep=""))

filter_1000genome_file = function(x, population){
  svp <- ScanVcfParam(info=c(population, "VT"), geno=NA) # only read in AF and VT
  vcf1 = readVcf(x, "hg19", svp)
  cat("readin",x)
  vcf2 = vcf1[sapply(info(vcf1)[[population]], max) > 0.01] # EUR_AF > 0.01, European pulation, there could be multiple minor allele, choose the largest one
  vcf2 = vcf2[isSNV(vcf2)]
  r = rowRanges(vcf2)
  cat("filter",x)
  return(r)
}

l         = lapply(IPF, filter_1000genome_file, population = "EUR_AF")
l         = do.call("c", l)



# get 1000 genome snps distance to tss
seqlevelsStyle(l) <- "UCSC"
dist      = distanceToNearest(l, tss_19)
l$dist    = mcols(dist)[,1]

# delete ARB.SNPS in 1000_genome_snps
f = findOverlaps(ARB.SNPS.hg19, l) # same number
l = l[-subjectHits(f)] # delete ARB SNPs in control SNPs


saveRDS(l, file = "~/project/variant-scores_manuscript/data/1KG_processed_eur.RDS")


KG_processed_eur <- readRDS("~/project/variant-scores_manuscript/data/1KG_processed_eur.RDS")
noncoding = readRDS("~/project/variant-scores_manuscript/data/noncoding_hg19.RDS")

f = findOverlaps(KG_processed_eur, noncoding)
a = KG_processed_eur[queryHits(f)]
saveRDS(a, file = "~/project/variant-scores_manuscript/data/1KG_processed_eur_noncoding.RDS")

a2 = readRDS("~/project/variant-scores-2/data/genetic/1KG_processed_eur_noncoding.RDS")
