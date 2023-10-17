library(GenomicFeatures)
library(GenomicRanges)
library(pegas)
library(biomaRt)
library(VariantAnnotation)
library(rtracklayer)
# Rcpp need >= 1.0.7


GENCODE_FILE.GR37 = "~/project/variant-scores_manuscript/external_data/gencode/rel36/gencode.v36lift37.annotation.gff3.gz"
db          = makeTxDbFromGFF(GENCODE_FILE.GR37, format=c("gff3"))
# exons       = exons(db, columns=c("tx_id", "tx_name"))
cds         = cds(db, columns=c("tx_id", "tx_name"))
# strand(exons) = "*"
strand(cds) = "*"
noncoding = gaps(cds)


saveRDS(noncoding, "~/project/variant-scores_manuscript/data/noncoding_hg19.RDS")

# noncoding2 = readRDS("~/project/variant-scores-2/data/genetic/noncoding_hg19.RDS")
# all.equal(noncoding, noncoding2)
# str(noncoding)
# str(noncoding2)