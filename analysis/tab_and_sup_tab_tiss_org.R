library(readr)

p_value = read_csv("./sup_data/sup_data_pairwise-tis-vs-org-individual.csv.gz")

p_value_pr = p_value[p_value$measure == "pr",]
p_value_pr = split( p_value_pr , f = p_value_pr$disease )
p_value_roc = p_value[p_value$measure == "roc",]
p_value_roc = split( p_value_roc , f = p_value_roc$disease )


get_stats = function(x){
  p = do.call("rbind", x)
  p_sig = p[p$p_value <0.05,]
  p_sig[p_sig$higher_median != "DHS-Weighted",]$higher_median = "org_score"
  ta = tapply(p_sig$higher_median, p_sig$score1, table)
  ta = do.call("rbind", ta)
  ta = as.data.frame(ta)
  ta$score = rownames(ta)
  ta$tie = 111-(ta$org_score + ta$`DHS-Weighted`)
  ta = ta[c(3,5,4,2,1),]
  ta1 = ta[,c("score", "DHS-Weighted", "org_score", "tie")]
  rownames(ta1) = NULL
  colnames(ta1)[2:3] = c("DHS_weighted_win", "DHS_weighted_lose")
  return(ta1)
}

o1_pr = get_stats(p_value_pr)
o2_roc = get_stats(p_value_roc)


write_csv(o1_pr, "./sup_data/tab_tiss_org.csv")
write_csv(o2_roc, "./sup_data/sup_tab_tiss_org.csv")
