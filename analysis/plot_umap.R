library(ggplot2)
library(ggrepel)


tt = readRDS("./data/analysis/coeffi_clusters.RDS")
tt$clst_nm = factor(tt$clst_nm, 
                    levels = c("cardiovascular disease/others", 
                               "digestive system disease/cancer",
                               "immune system disease",
                               "immune system disease/autoimmune disease",
                               "mental or behavioural disorder",
                               "heterogenous",
                               "skin cancer"))



# pick the representative diseases for each cluster
rep_ds = c("rheumatoid_arthritis", "autoimmune_disease", "inflammatory_bowel_disease", 
           "celiac_disease", "systemic_scleroderma", "acute_lymphoblastic_leukemia", 
           "cutaneous_melanoma", "non-melanoma_skin_carcinoma",
           "pancreatic_carcinoma", "diabetes_mellitus", "digestive_system_disease",
           "heart_failure", "alzheimer's_disease", "glaucoma", 
           "alcohol_dependence", "lung_carcinoma", "obesity", 
           "schizophrenia", "anxiety_disorder"
)

tt$disease2 = ""
tt$disease2[tt$disease %in% rep_ds] = tt$disease[tt$disease %in% rep_ds]

# generating 20 colors using https://medialab.github.io/iwanthue/
colors_20 = c("#30cb52",
              "#f76ff1",
              "#7fa900",
              "#9f7cff",
              "#e4c40a",
              "#623ba9",
              "#db5702",
              "#03b5ed",
              "#ff4797",
              "#00d7c1",
              "#f04ac4",
              "#007132",
              "#a1007d",
              "#e3c18d",
              "#ff97ec",
              "#6f5000",
              "#83c7f3",
              "#8c3836",
              "#be84aa",
              "#7e3b67")

color_cluster = setNames(colors_20[1:7], c(sort(unique(tt$clst_nm), decreasing = F))) # i don't want this to change
tt$color = color_cluster[match(tt$clst_nm, names(color_cluster))]


p <- ggplot(tt, aes(d1, d2, label = disease2, color = clst_nm)) +
  geom_point(size=2) +
  scale_x_continuous(expand=expansion(mult=1/15)) + scale_y_continuous(expand=expansion(mult=1/10)) +
  geom_text_repel(force=18,alpha=20,box.padding = 1.2, max.overlaps = Inf) +
  #cscale_color_brewer(palette="Dark2",name="cluster") +
  scale_color_manual(values=color_cluster, name = "cluster") +
  theme_minimal() + 
  theme(
    legend.position = c(1, 1),
    legend.justification = c("right", "top")
  ) + 
  labs(x = "UMAP1", y = "UMAP2")
cowplot::save_plot(filename="./plot/umap_subnames.pdf",
                   p, base_height = 6.5, base_width = 9)

p <- ggplot(tt, aes(d1, d2, label = disease, color=clst2_nm)) +
  geom_point(size=2) +
  # scale_color_brewer(palette="Dark2",name="cluster") +
  # scale_fill_manual(values=c("grey", "black"))+ 
  scale_color_grey() + 
  # theme_minimal() + 
  theme_classic()+ 
  labs(x = "UMAP1", y = "UMAP2", 
       color = "cluster")+
  theme(plot.title = element_text(size = 13, face = "bold"),
        legend.title=element_text(size=13), 
        legend.text=element_text(size=13))

p

cowplot::save_plot(filename="./plot/umap_two_clusters.pdf",p, base_height = 3, base_width = 4)
# plot group
# p <- ggplot(tt, aes(d1, d2, label = disease, color=clst_nm)) +
#   geom_point(size=3) +
#   scale_x_continuous(expand=expansion(mult=1/15)) + scale_y_continuous(expand=expansion(mult=1/10)) +
#   geom_text_repel(force=13,alpha=.75,box.padding = 2) +
#     scale_color_brewer(palette="Dark2",name="cluster") +
#   theme_minimal()

