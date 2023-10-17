library(ggplot2)

org_auc <- readRDS("./data/analysis/perf_org.RDS")

org_auc$match = factor(org_auc$match, levels = c("SNPsnap_TSS", "SNPsnap", "TSS", "random"))

# level_order <- c("random","SNPsnap_TSS", "SNPsnap", "TSS")
org_roc = org_auc[org_auc$measure == "roc",]
org_pr  = org_auc[org_auc$measure == "pr",]


# plot the performance for every disease, just look at them
p = ggplot(org_roc, aes(factor(match),auc,col=variant_score)) + geom_quasirandom() + 
  facet_wrap("PTPE") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
p = ggplot(org_pr, aes(factor(match),auc,col=variant_score)) + geom_quasirandom() + 
  facet_wrap("PTPE") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# plot normalize
nor_org = org_auc[org_auc$match != "random",]
nor_org$normalizeauc = nor_org$auc/(org_auc[org_auc$match == "random", "auc"])

nor_org_roc = nor_org[nor_org$measure == "roc",]
nor_org_pr  = nor_org[nor_org$measure == "pr",]

pdf(file="./plot/matching_strategy.pdf", width = 9, height = 3.5)
p = ggplot(nor_org_roc, aes(factor(match),normalizeauc, col = match)) + geom_quasirandom(size = .75) + 
  facet_wrap("variant_score", nrow = 1)
p + geom_hline(yintercept=1, linetype="dashed", 
               color = "red", size=0.3)+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  xlab("matching strategy") + 
  ylab("AUROC (normalize for random)") + 
  ggtitle("Performance of different matching strategy: normalize by random matching (AUROC)")




p = ggplot(nor_org_pr, aes(factor(match),normalizeauc, col = match)) + geom_quasirandom(size = .75) + 
  facet_wrap("variant_score", nrow = 1)
p + geom_hline(yintercept=1, linetype="dashed", 
               color = "red", size=0.3)+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  xlab("matching strategy") + 
  ylab("AUPRC (normalize for random)") + 
  ggtitle("Performance of different matching strategy: normalize by random matching (AUPRC)")
dev.off()
