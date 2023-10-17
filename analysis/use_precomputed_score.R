library(GenomicRanges)
library(IRanges)
library(rtracklayer)

# Step 1: download precompuated scores from: https://dataverse.harvard.edu/privateurl.xhtml?token=ead542d0-7527-4c55-bc34-bc50fc1c4cd5

# Step 2: have your input file ready as a txt file. 
## Format: chromosome:position. (genome assembly: GRCh37 - hg19)
## Please see example input files as: "examples.txt"
## we only have the scores for 1-22 chromosomes.

# Step 3: change the following into the path of the precomputed score and your input file

yoursample = "./data/examples.txt"
precomputed_score = "./precomputed_score_combined/"

# Step 4: prepare your sample to read in score
prepare_gr = function(x){
  sample = read.table(yoursample) # the location of your input file
  # get grange object
  chr = gsub("([0-9]+):([0-9]+)",
             "\\1",
             sample[,1])
  pos = as.integer(gsub("([0-9]+):([0-9]+)",
                        "\\2",
                        sample[,1]))
  gr2 <- GRanges(seqnames=chr, IRanges(pos,
                                       pos))
  seqlevelsStyle(gr2) <- "ucsc"
  return(gr2)
}

sample_gr = prepare_gr(yoursample)

# Step 5: if you want to read in file for one specific disease, do the following.
## remember to specify the disease name first
diseasename = "acute_lymphoblastic_leukemia" # change disease name into your interested ones
gr_disease = import(paste0(precomputed_score, diseasename, ".bw"), # the location of the downloaded file
                    which = sample_gr)

# Step 6: if you want to read in all scores: do the following
read_all_scores = function(file_path, gr2){
  files = Sys.glob(paste0(precomputed_score, "*.bw"))
  gr = sapply(files, function(x){
    name = gsub(".bw","",basename(x))
    print(name)
    gr = import(x, which = gr2)
    names(mcols(gr)) = name
    return(mcols(gr)[name])
  })
  names(gr) = NULL
  gr = do.call("cbind", gr)
  mcols(gr2) = gr
  return(gr2)
}

gr_all_diseases = read_all_scores(precomputed_score, sample_gr)


