## function to get control snps using bins in control snps
get_control_snps_controlbins = function(Disease_SNPS, control_SNPS, bin = 50, times = 10){
  set.seed(0525)
  q = quantile(control_SNPS$dist, probs = seq(0, 1, length.out = bin+1), na.rm = FALSE,
               names = TRUE, type = 7)
  control.snps = vector("list", bin)
  for (i in 1:bin){
    d = Disease_SNPS[Disease_SNPS$dist < q[i+1] & Disease_SNPS$dist >= q[i]]
    if (length(d) != 0){
      ctrl              = sample(control_SNPS[control_SNPS$dist < q[i+1] & control_SNPS$dist >= q[i]], length(d)*10)
      ctrl2             = GRanges(granges(ctrl), 
                                  type = 0, 
                                  snpID = paste(gsub("chr","",seqnames(ctrl)), start(ctrl), sep = ":"), 
                                  rsID = names(ctrl), 
                                  input_snp = rep(d$input_snp, each = 10), 
                                  dist = ctrl$dist)
      control.snps[[i]] = ctrl2
      
    }else{
      # control.snps[[i]] = control_SNPS[Disease_SNPS$dist < q[i+1] & Disease_SNPS$dist >= q[i]] # this is length zero
      control.snps[[i]] = granges(control_SNPS[Disease_SNPS$dist < q[i+1] & Disease_SNPS$dist >= q[i]])
    }
    print(i)
  }
  control.snps = do.call("c", control.snps)
  return(control.snps)
}
