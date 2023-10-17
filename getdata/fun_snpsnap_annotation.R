library(reshape2)

# read in ld0.8 file download from snpsnap database
ld0.8 = read.table(gzfile("~/project/variant-scores_manuscript/snpsnap/ld0.8_collection.tab.gz"), sep = "\t", header = T)

snpsnap_annotation = function(output_file){
  matched_file = read.table(output_file, sep = "\t", header = T, stringsAsFactors = F)
  
  ### they already exclude it when running on SNPsnap website
  # exclude HLA regions 6:25000000-6:35000000
  # https://data.broadinstitute.org/mpg/snpsnap/match_snps.html
  # HLA regions explicitly mentioned on the website above: 
  # Go to 'SNP exclusions', point to the tick before 'exclude HLA SNPs', the exclusion criteria will show up
  split.match = strsplit(matched_file$lead_snp, ":")
  chr = sapply(split.match, function(x) x[1])
  pos = as.numeric(sapply(split.match, function(x) x[2]))
  exclud = chr == 6 & pos >=25000000 & pos <= 35000000
  sort(matched_file[exclud,]$snps_to_match) # look at if it is correct
  matched_file_filtered = matched_file[!exclud,]
  
  
  # from wide format to long format
  matched = reshape2::melt(matched_file_filtered, id.vars=c("lead_snp"))
  matched$variable = gsub("Set_", "", matched$variable) 
  colnames(matched) = c("input_snp", "set", "snpID") # set colnames to be the same with SNPsnap website
  matched = matched[order(matched$input_snp),]
  matched = matched[,c("set", "input_snp", "snpID")] 
  
  
  # match for annotation
  input = data.frame(snpID = matched[matched$set == 1,"input_snp"])
  mi = match(input$snpID, ld0.8$snpID)
  mo = match(matched$snpID, ld0.8$snpID)
  input_anno = ld0.8[mi,]
  output_anno = ld0.8[mo,]
  output_anno = cbind(matched[,c("set", "input_snp")], output_anno)
  return(list(input_anno, output_anno))
}
