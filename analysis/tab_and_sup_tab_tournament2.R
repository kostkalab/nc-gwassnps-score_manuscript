library(readr)

p_value = read_csv("./sup_data/sup_data_pairwise-tis-weighted-aggregated.csv.gz")
p_value_by_disease = read_csv("./sup_data/sup_data_pairwise-tis-weighted-individual.csv.gz")

p_value$lower_median = apply(p_value, 1, function(y){
  return(setdiff(c(y[1], y[3]), y[7]))
})

p_value$higher_median = factor(p_value$higher_median, levels = c("DHS-Weighted", "Genoskyline-Weighted", "Fitcons2-Weighted"))
p_value$lower_median = factor(p_value$lower_median, levels = c("DHS-Weighted", "Genoskyline-Weighted", "Fitcons2-Weighted"))

p_value_sig = p_value[p_value$p_value < 0.05,]


t1 = tapply(p_value_sig$higher_median, p_value_sig$measure, table)
t1 = t(do.call("rbind", t1)) # wins
t2 = tapply(p_value_sig$lower_median, p_value_sig$measure, table)
t2 = t(do.call("rbind", t2)) # loses


p_value_by_disease$lower_median = apply(p_value_by_disease, 1, function(y){
  return(setdiff(c(y[1], y[3]), y[8]))
})
p_value_by_disease$higher_median = factor(p_value_by_disease$higher_median, levels = c("DHS-Weighted", "Genoskyline-Weighted", "Fitcons2-Weighted"))
p_value_by_disease$lower_median = factor(p_value_by_disease$lower_median, levels = c("DHS-Weighted", "Genoskyline-Weighted", "Fitcons2-Weighted" ))
p_value_by_disease_sig = p_value_by_disease[p_value_by_disease$p_value < 0.05,]


t3 = tapply(p_value_by_disease_sig$higher_median, p_value_by_disease_sig$measure, table)
t3 = t(do.call("rbind", t3)) # wins
t4 = tapply(p_value_by_disease_sig$lower_median, p_value_by_disease_sig$measure, table)
t4 = t(do.call("rbind", t4)) # loses

format_df = function(df_pr){
  df_pr = data.frame(df_pr)
  df_pr[,5] = apply(df_pr, 1, function(x) 2-x[1]-x[2])
  df_pr[,6] = apply(df_pr, 1, function(x) 222-x[3]-x[4])
  df_pr$score = rownames(df_pr)
  df_pr = df_pr[,c(7,1,2,5,3,4,6)]
  colnames(df_pr) = c("score" ,"win-aggre", "lose-aggre", "tie-aggre", "win", "lose", "tie")
  df_pr = df_pr[order(df_pr[,5], decreasing = T),]
  rownames(df_pr) = NULL
  return(df_pr)
}

df_pr = cbind(t1[,1],t2[,1],t3[,1], t4[,1])
df_pr = format_df(df_pr)
df_roc = cbind(t1[,2], t2[,2], t3[,2], t4[,2])
df_roc = format_df(df_roc)


write_csv(df_pr, "./sup_data/tab_tournament2.csv")
write_csv(df_roc, "./sup_data/sup_tab_tournament2.csv")
