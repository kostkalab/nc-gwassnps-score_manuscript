library(ComplexHeatmap)

# get coefficients and epigenome annotation
#------------------------------------------

coeff              <- readRDS("./data/analysis/coefficients.RDS")
top_mean           <- sapply(coeff, function(x) {x[,1]})
rownames(top_mean) <- rownames(coeff[[1]])
coeffi             <-  top_mean[-c(1,2),] # remove rowmeans and intercept, only concentrate on tissues

# normalize coeffi
coeffi_norm = apply(coeffi,2,function(x) {x=x-min(x); x=  x/quantile(x,.95); x[x>1]=1; return(x)}) # QQ: kind of normalize it; anything larger than 1 is 1
co_scale = coeffi_norm

## read in cluster name, and top tissue file
tt = readRDS("./data/analysis/coeffi_clusters.RDS")
clusterdf = readRDS("./data/analysis/coeffi_top5tissues.RDS")
clusterdf$cluster_nm_short = tt$clst_shortnm[match(clusterdf$cluster_id, tt$clst)] # add short name to plot, to prevent overlapping


# retrieve ANATOMY, TYPE, and GROUP for 127 epigenomes. 
library(erma) # https://bioconductor.org/packages/devel/bioc/vignettes/erma/inst/doc/erma.html
meta          <- mapmeta()
met           <- as.data.frame(meta[,c("Epigenome.ID..EID.", "GROUP","Standardized.Epigenome.name", "ANATOMY", "TYPE")])
rownames(met) <- met[,1]
met           <- met[rownames(coeffi),]




library(ComplexHeatmap)
library(circlize)


dfhp = co_scale[clusterdf$tissue,] # get the table for the top tissue heatmap, select the tissues in the top tissues

## prepare to color the group
# get a vector that stores disease and the cluster it belongs to
sub_grp = tt$clst_shortnm
names(sub_grp) = tt$disease

# get a vector that stores the tissue and the cluter it belongs to
sub_grp3 = clusterdf$cluster_nm_short
names(sub_grp3) = clusterdf$tissue

m = match(names(sub_grp3), met$Epigenome.ID..EID.) # match to add the real tissue name (group name) instead of the tissue id


# define color for heatmap
col_fun = colorRamp2(c(0, 0.5, 0.75,1), c("white", "#F6D55C", "#EB3C13", "#AB0000"))

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

pie(rep(1, 20), col=colors_20) # look at those colors

color_cluster = setNames(colors_20[1:7], c(sort(unique(tt$clst_nm)))) # define a color to each cluster
color_tiss    = setNames(colors_20[c(8:9,5,10:19,7,20)], c(sort(unique(met$ANATOMY[m])))) # define a color to each tissue group

## make the annotation object
hleft = rowAnnotation(clus_nm = tt$clst_nm, 
                      col = list(clus_nm = color_cluster),
                      show_legend = c("tiss_nm" = FALSE), # don't show legend
                      show_annotation_name = FALSE
)

hbottom = columnAnnotation(tiss_nm = met$ANATOMY[m], 
                           col = list(tiss_nm = color_tiss), 
                           show_legend = c("tiss_nm" = FALSE), # don't show legend,
                           show_annotation_name = FALSE
)




col_labels = structure(paste(met$Epigenome.ID..EID., met$ANATOMY, sep = "-"), names = met$Epigenome.ID..EID.)

# read in the performance data to add the performance annotation
DHS_logreg = readRDS("./data/analysis/perf_DHS_logreg.RDS")
DHS_logreg = DHS_logreg[DHS_logreg$measure == "pr",] # only keep PR
DHS_logreg = data.frame(t(tapply(DHS_logreg$auc, list(DHS_logreg$regularization, DHS_logreg$PTPE), median)))
DHS_logreg$diff = DHS_logreg$min - DHS_logreg$baseline

m2 = match(colnames(dfhp), rownames(DHS_logreg))
DHS_logreg = DHS_logreg[m2,] # organize so the sorting is the same as the heatmap.


col_fun_r = colorRamp2(c(0.09, 0.36), c( "white", "blue"))
row_right = rowAnnotation(AUPRC = DHS_logreg$min
                          # , diff = DHS_logreg$diff
                          , col = list(AUPRC = col_fun_r))


pdf(file="./plot/heatmap_toptiss.pdf", width = 8, height = 6)
hp = Heatmap(t(dfhp), 
             col = col_fun, 
             name = "coefficients",
             row_split = sub_grp, column_split = sub_grp3,
             cluster_columns = F, 
             show_row_names = F,
             column_names_gp = gpar(fontsize = 8), 
             row_title_gp = gpar(fontsize = 8), 
             column_title_gp = gpar(fontsize = 8), 
             row_title_rot = 0,
             column_labels = col_labels[colnames(t(dfhp))],
             left_annotation = hleft, 
             cluster_rows = F, 
             bottom_annotation = hbottom, 
             right_annotation = row_right
)
draw(hp, column_title = "cluster-specific tissues")
dev.off()


## draw the one with row names on it as supplemental figure
pdf(file="./plot/heatmap-toptiss-rownames.pdf", width = 10, height = 14)
hp = Heatmap(t(dfhp), 
             col = col_fun, 
             name = "coefficients",
             row_split = sub_grp, column_split = sub_grp3,
             cluster_columns = F, 
             show_row_names = T,
             column_names_gp = gpar(fontsize = 8), 
             row_names_gp = gpar(fontsize = 8),
             row_title_gp = gpar(fontsize = 8), 
             column_title_gp = gpar(fontsize = 8), 
             row_title_rot = 0,
             column_labels = col_labels[colnames(t(dfhp))],
             left_annotation = hleft, 
             cluster_rows = F,
             bottom_annotation = hbottom
)
draw(hp, column_title = "cluster-specific tissues")
dev.off()

## heatmap_data = list(top_matrix = dfhp, 
##                     toptissud = sub_grp3, 
##                     cluster = sub_grp)
## 
## saveRDS(heatmap_data, "~/project/variant-scores-2/analysis/GWAS_catalog_3.0/aggre/result_data/heatmap_data.RDS")
## 
