library(GenomicFeatures)
library(GenomicRanges)
library(pegas)
library(biomaRt)
library(VariantAnnotation)
library(rtracklayer)

FILE_CHAIN   ="~/project/variant-scores_manuscript/external_data/ucsc/hg38ToHg19.over.chain"
# GENCODE_FILE.GR38 = "/data/projects/annotation/GENCODE/rel29/gencode.v29.annotation.gff3.gz"
GENCODE_FILE.GR38 = "~/project/variant-scores_manuscript/external_data/gencode/rel36/gencode.v36.annotation.gff3.gz"


#------------ 
# GET TSS FOR PROTEIN CODING GENES
# filter for protein coding transcripts
db          = makeTxDbFromGFF(GENCODE_FILE.GR38, format=c("gff3"))
transcripts = transcripts(db, columns=c("tx_id", "tx_name"))

ensembl             = useMart("ensembl",dataset="hsapiens_gene_ensembl") # hg38
out                 = getBM(attributes=c("ensembl_transcript_id_version", "transcript_biotype"), 
                            filters="ensembl_transcript_id_version", values=transcripts$tx_name, mart=ensembl) #this runs for a very long time
out                 = out[match(transcripts$tx_name, out$ensembl_transcript_id_version),]

transcripts$biotype = out$transcript_biotype
transcripts         = transcripts[!is.na(transcripts$biotype)]
transcripts_protein = transcripts[transcripts$biotype == "protein_coding"]
tss_38              = resize(transcripts_protein, width=1, fix='start')
chain               = import.chain(FILE_CHAIN)
tss_19              = unlist(GRangesList(liftOver(tss_38,chain)))


saveRDS(tss_38, "~/project/variant-scores_manuscript/data/tss_38.RDS")
saveRDS(tss_19, "~/project/variant-scores_manuscript/data/tss_19.RDS")