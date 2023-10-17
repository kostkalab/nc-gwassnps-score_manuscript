library(ggplot2)
library(data.table)
library(ggpubr)
# plot DIVAN, DHS-logreg, and GWAVA, not plot DHS-baseline


# organism score
org_score = "GenoCanyon"

# auc_all = readRDS("~/projects/variant-scores-2/data/DIVAN/data/auc_all_GWAVA.RDS")
auc_all = readRDS("./data/analysis/perf_all_divan.RDS")
# exclude body weight from auc_all
exclude = "body_weight"
auc_all = auc_all[!auc_all$PTPE %in% exclude,]

auc = tapply(auc_all$auc, list(auc_all$measure, auc_all$PTPE, auc_all$type), median)
auc_pr = data.frame(t(auc[,,"pr"]))

p1 <- ggplot(auc_pr, aes(x=GenoCanyon, y=DIVAN)) +
  geom_point() + 
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
  labs(x = "GenoCanyon", y = "DIVAN")
#xlim(0.085,0.24) + 
#ylim(0.09,0.27)

p2 <- ggplot(auc_pr, aes(x=GenoCanyon, y=min)) +
  geom_point() + 
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
  labs(x = "GenoCanyon", y = "DHS-weighted") 
#xlim(0.085,0.24) + 
#ylim(0.09,0.27)

p3 <- ggplot(auc_pr, aes(x=DIVAN, y=min)) +
  geom_point() + 
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
  labs(x = "DIVAN", y = "DHS-weighted")
#xlim(0.085,0.24) + 
#ylim(0.09,0.27)


pdf(file="./plot/perf_divan_scatter_plot.pdf", width = 9, height = 3)
ggarrange(p1, p2, p3, 
          # labels = c("A", "B", "C"),
          nrow = 1)+
  theme(plot.margin = margin(0.1,0.3,0.1,0.1, "cm")) 
dev.off()
