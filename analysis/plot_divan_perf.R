library(ggplot2)
library(data.table)
library(dplyr)
# plot DIVAN, DHS-logreg, and GWAVA, not plot DHS-baseline


# organism score
org_score = "GenoCanyon"

# auc_all = readRDS("~/projects/variant-scores-2/data/DIVAN/data/auc_all_GWAVA.RDS")
auc_all = readRDS("./data/analysis/perf_all_divan.RDS")
p_value = readRDS("./data/analysis/pairwise_DIVAN_top.RDS")

# p_value = do.call("rbind", p_value)
# get the average auc
auc_mean = reshape2::melt(tapply(auc_all$auc, list(auc_all$type, auc_all$PTPE, auc_all$measure), mean))
colnames(auc_mean) = c("type", "PTPE", "measure", "mean")
auc_min  = reshape2::melt(tapply(auc_all$auc, list(auc_all$type, auc_all$PTPE, auc_all$measure), min))
colnames(auc_min) = c("type", "PTPE", "measure", "min")
auc_max  = reshape2::melt(tapply(auc_all$auc, list(auc_all$type, auc_all$PTPE, auc_all$measure), max))
colnames(auc_max) = c("type", "PTPE", "measure", "max")
auc_df = auc_mean %>% full_join(auc_max, by = c("type", "PTPE", "measure")) %>% 
  full_join(auc_min, by = c("type", "PTPE", "measure"))


# Default bar plot
# roc

# rank
a = auc_df[auc_df$type == "pr" & auc_df$measure %in% c("min", "DIVAN", org_score),]
firstup = function(x){
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
a$PTPE = firstup(gsub("_", " ",a$PTPE))

# have a column states that org, DIVAN and our which is better
# three categories
# linsight best, DIVAN best, our best
a = data.table(a)

a[measure == "min",      measure := "Tissue-Weighted (DHS)"]


df_comp = a[,c("PTPE", "measure", "mean")]
df_comp_r = reshape(df_comp, varying = list(c("DIVAN", org_score, "Tissue-Weighted (DHS)")), 
                    v.names = "mean", 
                    timevar = "measure", 
                    idvar = "PTPE",
                    direction = "wide")
id = apply(df_comp_r[,2:4], 1,function(x){
  m = c("DIVAN", org_score, "Tissue-Weighted (DHS)")[which.max(x)]
  return(m)
})
a = as.data.frame(a)
a$id = id


# plot




# rank
a$measure = droplevels(a$measure)
a_best = a[a$measure == a$id, ]

rank = a_best$PTPE[order(a_best$mean)]
a$PTPE = factor(a$PTPE, levels = rank)

# pick the diseaes that have a p value larger than 0.05
disease = p_value[p_value$p_value >= 0.05,"disease"]
a$id = as.character(a$id)
a[a$PTPE %in% disease, "id"] = "non-sig"
a$id = factor(a$id, levels = c(org_score, "DIVAN", "Tissue-Weighted (DHS)", "non-sig"),
              labels = c(org_score, "DIVAN", "DHS Tissue-weighted", "non\n-sig"))
a$measure = factor(a$measure, levels = c(org_score, "DIVAN", "Tissue-Weighted (DHS)"),
                   labels = c(org_score, "DIVAN", "DHS Tissue-weighted"))




pdf("./plot/perf_divan.pdf", height = 5, width = 13)
p<- ggplot(a, aes(x=PTPE, y=mean, fill=measure)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=min, ymax=max), width=.2,
                position=position_dodge(.9))  + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95)) + 
  facet_grid(. ~ id, scales = "free_x", space = "free_x") + 
  geom_hline(yintercept=1/11, linetype="dashed") + 
  theme(# The new stuff
    strip.text = element_text(size = 13), 
    axis.text = element_text( size = 13 ),
    axis.text.x = element_text( size =11 ),
    axis.title = element_text( size = 13 ),
    legend.title=element_text(size=14), 
    legend.text=element_text(size=13),
    strip.background = element_blank(),
    plot.margin = margin(b=5,r=5,t=7,l=5, "pt")) + 
  labs(x = NULL, y = "Area under the PR Curve", fill = "Variant score")
print(p)
dev.off()



