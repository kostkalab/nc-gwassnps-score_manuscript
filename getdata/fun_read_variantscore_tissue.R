library(GenomicRanges)

# gr.hg19 is a GRange object with rownames rsID number
read_genoskyline = function(gr.hg19,
                            DIR_GENO="./external_data/variant_scores/genoSkyline/BedGraph/Roadmap_E*.bedGraph.gz"){
  # this works well if there are duplicate snps; checked Jan 2020
  files = Sys.glob(DIR_GENO)
  basename(files)
  names = substr(basename(files),start=9,stop=12)
  drop  = which(names %in% c("E035","E123"))
  names = names[-drop]
  files = files[-drop]
  
  smat  = matrix(NA,ncol=length(names),nrow=length(gr.hg19))
  colnames(smat) = names
  rownames(smat) = names(gr.hg19)
  
  for(i in 1:length(names)){
    tmp = import(files[i],format="bedGraph") # no overlaping granges
    ov  = findOverlaps(gr.hg19,tmp)
    smat[queryHits(ov),colnames(smat) == names[i]] = tmp$score[subjectHits(ov)]
    print(names[i])
  }
  smat = as.data.frame(cbind(smat, mcols(gr.hg19)))
  return(smat)
}


read_DHS = function(snpset.locs,
                    DIR_DHS = "./external_data/avocado/DHS/data/avocado.Dnase.bw/*.avocado.bw"){
  snpset.locs     = snpset.locs[!(seqnames(snpset.locs) == "chrX")]
  files = Sys.glob(DIR_DHS)
  names = gsub(".DNase.allchrs.avocado.bw", "", basename(files))
  smat  = matrix(NA,ncol=length(names),nrow=length(snpset.locs))
  colnames(smat) = names
  rownames(smat) = names(snpset.locs)
  
  for(i in 1:length(names)){
    tmp = import(files[i], which = disjoin(snpset.locs))
    ov  = findOverlaps(snpset.locs,tmp)
    smat[queryHits(ov),colnames(smat) == names[i]] = tmp$score[subjectHits(ov)]
    print(names[i])
  }
  smat = as.data.frame(cbind(smat, mcols(snpset.locs)))
  return(smat)
}



read_Fitcons2 = function(snpset.locs,
                         DIR_FITCONS2 = "./external_data/variant_scores/fitcons2/data/*sco.bw"){
  files = Sys.glob(DIR_FITCONS2)
  names = gsub("-sco.bw", "", basename(files))
  smat  = matrix(NA,ncol=length(names),nrow=length(snpset.locs))
  colnames(smat) = names
  rownames(smat) = names(snpset.locs)
  
  for(i in 1:length(names)){
    tmp = import(files[i], which = disjoin(snpset.locs))
    ov  = findOverlaps(snpset.locs,tmp)
    smat[queryHits(ov),colnames(smat) == names[i]] = tmp$score[subjectHits(ov)]
    print(names[i])
  }
  smat = as.data.frame(cbind(smat, mcols(snpset.locs)))
  return(smat)
}


