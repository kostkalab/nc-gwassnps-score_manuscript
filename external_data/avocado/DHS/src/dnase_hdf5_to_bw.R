## make .hdf5 avocado Dnase file into .bw file

library(rhdf5)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(stringr) 

# get all file names
IPD = "/data/projects/annotation/avocado/data/"
IPF = Sys.glob(paste(IPD,"*.npz.hdf5",sep=""))

# get tissue and chromosome names
a1 = do.call(rbind, strsplit(gsub(".avocado.npz.hdf5|.DNase", "", basename(IPF)),
                             split="\\."))
tissues = unique(a1[,1])
chrs = unique(a1[,2])

# make them into granges files, concatenate them by tissue, 
# and write them out seperately by tissue
tmp = lapply(tissues, function(tiss){ # looping over tissue
  this.files = IPF[grep(str_c(tiss, "\\.DNase"),  IPF)] # get all file names one tissue
  this.RESA = apply(as.matrix(this.files), 1, function(x) h5read(x, "arr_0")) # read in all files in a list for one tissue
  ## is a list, of length 22; now set the names of list
  names(this.RESA) = gsub(".avocado.npz.hdf5|E[0-9]+.DNase.", "", basename(this.files))
  # make them into grange objects
  out = sapply(seq_along(this.RESA), function(ind){
    this.chr = names(this.RESA)[ind]
    gr = GRanges(seqnames = this.chr, 
                 ranges = IRanges(start = seq(1, as.integer(seqlengths(Hsapiens)[this.chr]/25*25-24), by = 25),
                                  end = seq(25, as.integer(seqlengths(Hsapiens)[this.chr]/25*25), by = 25)),
                 strand = "*", 
                 score = as.vector(this.RESA[[ind]]))
    return(gr)
  })
  towrite = do.call(c, out) # concatenate
  seqlengths(towrite)[sort(names(seqlengths(Hsapiens)[1:22]))] = 
    seqlengths(Hsapiens)[1:22][sort(names(seqlengths(Hsapiens)[1:22]))] # add seqlengths information
  export.bw(towrite, str_c("~/variant-scores/data/avocado.DNase/", tiss, 
                           ".DNase.allchrs.avocado.bw")) # export files
})
