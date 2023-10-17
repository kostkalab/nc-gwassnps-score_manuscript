library(readr)

p_value_by_disease = read_csv("./sup_data/sup_data_pairwise-divan-individual.csv.gz")
p_value_by_disease = read_csv("/sup_data_pairwise-divan-individual.csv.gz")
n = nrow(p_value_by_disease)/6 # how many diseases

p_value_by_disease$lower_median = apply(p_value_by_disease, 1, function(y){
  return(setdiff(c(y[1], y[3]), y[8]))
})

p_value_by_disease_sig = p_value_by_disease[p_value_by_disease$p_value < 0.05,]
t3 = tapply(p_value_by_disease_sig$higher_median, p_value_by_disease_sig$measure, table)
t3 = t(do.call("rbind", t3)) # wins
t4 = tapply(p_value_by_disease_sig$lower_median, p_value_by_disease_sig$measure, table)
t4 = t(do.call("rbind", t4)) # loses


format_df = function(df_pr){
  df_pr = data.frame(df_pr)
  df_pr[,3] = apply(df_pr, 1, function(x) 2*n-x[1]-x[2])
  df_pr$score = rownames(df_pr)
  df_pr = df_pr[,c(4,1,2,3)]
  colnames(df_pr) = c("score" ,"win", "lose", "tie")
  df_pr = df_pr[order(df_pr[,2], decreasing = T),]
  rownames(df_pr) = NULL
  return(df_pr)
}

df_pr = cbind(t3[,1], t4[,1])
df_pr = format_df(df_pr)
#df_pr$score = c("DHS-weighted", "GenoCanyon", "DIVAN")  ######### be really careful with this. 

df_roc = cbind(t3[,2], t4[,2])
df_roc = format_df(df_roc)
#df_roc$score = c("GenoCanyon", "DIVAN", "DHS_weighted")  ######### be really careful with this. 

write_csv(df_pr,  "~/project/variant-scores-2/data/GWAS_catalog_project_3.0/tables/tab_divan.csv")
write_csv(df_roc, "~/project/variant-scores-2/data/GWAS_catalog_project_3.0/tables/sup_tab_divan.csv")
