library(ggplot2)
library(tibble)
library(dplyr)
library(ggbeeswarm)
library(reshape2)
library(plyr)
library(data.table)


# read in SNPSNAP tissue specific scores
DHS_SNPSNAP_logreg <- readRDS("./data/analysis/perf_DHS_logreg.RDS")
Fit_SNPSNAP_logreg <- readRDS("./data/analysis/perf_Fit_logreg.RDS")
geno_SNPSNAP_logreg <- readRDS("./data/analysis/perf_geno_logreg.RDS")


# combine them together
all_p_SNPSNAP = rbind(DHS_SNPSNAP_logreg, Fit_SNPSNAP_logreg, geno_SNPSNAP_logreg)

# get the median of all 30 replicates
md = tapply(all_p_SNPSNAP$auc, 
            INDEX = list(all_p_SNPSNAP$measure, all_p_SNPSNAP$PTPE, 
                         all_p_SNPSNAP$variant_score, all_p_SNPSNAP$regularization), median)
md1 = reshape2::melt(md)
colnames(md1) = c("measure", "PTPE", "variant_score", "regularization", "auc")
md2 = md1[md1$regularization %in% c("min", "baseline"), ]
md2$regularization = factor(md2$regularization, levels = c("baseline", "min"))
md2$regularization = revalue(md2$regularization, c("baseline" = "baseline","min"="logreg"))


md2_pr = md2[md2$measure == "pr",]
md2_pr = data.table(md2_pr)
md2_pr[regularization == "logreg", regularization := "Tissue_Weighted"]
md2_pr[regularization == "baseline",      regularization := "Tissue_Mean"]
md2_pr$regularization = factor(md2_pr$regularization, levels = c("Tissue_Mean", "Tissue_Weighted"))

y1 = 0.36
y2 = 0.34
y3 = 0.32
y4 = 0.07
y5 = 0.06

p_dhs = 6.8e-20
p_fit = 3.7e-19
p_gno = 1.2e-19
p_dhs_fit = 5.3e-18
p_dhs_gno = 7.0e-11

pdf(file="./plot/perf_tiss.pdf", width = 5.5, height = 6)
p = ggplot(md2_pr,       aes(regularization, auc,col = variant_score)) + geom_boxplot() + 
  # ggtitle("Performance of Tissue_Mean and Tissue_Weighted") + 
  geom_quasirandom(dodge.width=0.75, size = 0.8) + 
  geom_hline(yintercept=1/11, linetype="dashed") + 
  ylab("Average precision") + 
  xlab(NULL)+
  theme_bw() + 
  theme(legend.position = "none", 
        # The new stuff
        strip.text = element_text(size = 13), 
        axis.text = element_text( size = 13 ),
        axis.text.x = element_text( size =13 ),
        axis.title = element_text( size = 13 ))+
  scale_x_discrete(labels=c("Tissue_Mean" = "Tissue-mean", 
                            "Tissue_Weighted" = "Tissue-weighted"))

p + annotate(geom="text",x=1.5,y=y1+0.01, label=paste0("p.val = ", p_fit))  + # fitcons
  annotate(geom="segment",x="Tissue_Mean",xend="Tissue_Weighted",y=y1,yend=y1) +
  annotate(geom="segment",x="Tissue_Mean",xend="Tissue_Mean",y=y1,yend=y1-0.01) +
  annotate(geom="segment",x="Tissue_Weighted",xend="Tissue_Weighted",y=y1,yend=y1-0.01) +
  annotate(geom="text",x=1.3,y=y2+0.01, label=paste0("p.val = ", p_dhs)) + # DHS
  annotate(geom="segment",x=0.75,xend=1.75,y=y2,yend=y2) +
  annotate(geom="segment",x=0.75,xend=0.75,y=y2,yend=y2-0.01) +
  annotate(geom="segment",x=1.75,xend=1.75,y=y2,yend=y2-0.01) + 
  annotate(geom="text",x=2,y=y3+0.01, label=paste0("p.val = ", p_gno)) + # genoskyline
  annotate(geom="segment",x=1.25,xend=2.25,y=y3,yend=y3) +
  annotate(geom="segment",x=1.25,xend=1.25,y=y3,yend=y3-0.01) +
  annotate(geom="segment",x=2.25,xend=2.25,y=y3,yend=y3-0.01) + 
  
  annotate(geom="text",x=2,y=y4+0.01, label=paste0("p.val = ", p_dhs_fit)) + # DHS vs Fitcons
  annotate(geom="segment",x=1.75,xend=2,y=y4,yend=y4) +
  annotate(geom="segment",x=1.75,xend=1.75,y=y4,yend=y4+0.005) +
  annotate(geom="segment",x=2,xend=2,y=y4,yend=y4+0.005) + 
  
  annotate(geom="text",x=2.25,y=y5-0.01, label=paste0("p.val = ", p_dhs_gno)) + # DHS vs geno
  annotate(geom="segment",x=1.75,xend=2.25,y=y5,yend=y5) +
  annotate(geom="segment",x=1.75,xend=1.75,y=y5,yend=y5+0.005) +
  annotate(geom="segment",x=2.25,xend=2.25,y=y5,yend=y5+0.005)
dev.off()


# plot scatter plot
md2_pr_wide = reshape(md2_pr, timevar = "regularization", idvar = c("PTPE", "variant_score", "measure"),direction = "wide")

pdf(file="./plot/perf_tiss_scatter.pdf", width = 6, height = 6)
p <- ggplot(md2_pr_wide, aes(x=auc.Tissue_Mean, y=auc.Tissue_Weighted, col = variant_score)) +
  geom_point() + 
  facet_wrap("variant_score", nrow = 2) + 
  theme_minimal()+ 
  theme(legend.position = c(0.75, 0.2), 
        # The new stuff
        strip.text = element_text(size = 13), 
        axis.text = element_text( size = 13 ),
        axis.text.x = element_text( size =14 ),
        axis.title = element_text( size = 13 ),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=13)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  labs(x = "Tissue-mean", y = "Tissue-weighted", col = "Variant score")
p
dev.off()
