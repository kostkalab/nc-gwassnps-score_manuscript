library(ComplexHeatmap)
library(circlize)
# get coefficients and epigenome annotation
#------------------------------------------

coeff              <- readRDS("./data/analysis/coefficients.RDS")
top_mean           <- sapply(coeff, function(x) {x[,1]})
rownames(top_mean) <- rownames(coeff[[1]])
coeffi             <-  top_mean[-c(1,2),] # remove rowmeans and intercept, only concentrate on tissues

# retrieve ANATOMY, TYPE, and GROUP for 127 epigenomes. 
library(erma) # https://bioconductor.org/packages/devel/bioc/vignettes/erma/inst/doc/erma.html
meta          <- mapmeta()
met           <- as.data.frame(meta[,c("Epigenome.ID..EID.", "GROUP","Standardized.Epigenome.name", "ANATOMY", "TYPE")])
rownames(met) <- met[,1]
met           <- met[rownames(coeffi),]

tt = readRDS("./data/analysis/coeffi_clusters.RDS")
clusterdf = readRDS("./data/analysis/coeffi_top5tissues.RDS")
# clusterdf$cluster_nm_short = tt$clst_shortnm[match(clusterdf$cluster_id, tt$clst)] # add short name to plot, to prevent overlapping



coeffi_norm = apply(coeffi,2,function(x) {x=x-min(x); x=  x/quantile(x,.95); x[x>1]=1; return(x)}) # QQ: kind of normalize it; anything larger than 1 is 1
# try z scores
# coeffi_norm = apply(coeffi,2,function(x) {x = (x-mean(x))/sd(x)})
co_scale = coeffi_norm

rs = apply(co_scale,1, function(rw) {
  rw2 = data.frame(rw=rw, group = tt$clst_shortnm) # add group name
  gm = tapply(rw2$rw, rw2$group, mean)# get the group mean
  associated_cluster = names(gm[which.max(gm)])
  return(associated_cluster)
})

# rename the five clusters using the short names


rs3 = rs

## get a vector that stores cluster name and disease
sub_grp = tt$clst_shortnm
names(sub_grp) = tt$disease

# col_fun = colorRamp2(c(0, 0.5, 1), c("white", "#F6D55C", "red"))
col_fun = colorRamp2(c(0, 0.5, 0.75,1), c("white", "#F6D55C", "#EB3C13", "#AB0000"))

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
color_cluster = setNames(colors_20[1:7], c(sort(unique(tt$clst_nm)))) # i don't want this to change
hleft = rowAnnotation(clus_nm = tt$clst_nm, 
                      col = list(clus_nm = color_cluster),
                      show_legend = c("tiss_nm" = FALSE), # don't show legend
                      show_annotation_name = FALSE
)

col_labels = structure(paste(met$Epigenome.ID..EID., met$ANATOMY, sep = "-"), names = met$Epigenome.ID..EID.)

pdf(file="./plot/heatmap_alltissues.pdf", width = 12, height = 7)
hp = Heatmap(t(co_scale), 
             col = col_fun, 
             name = "scaled coefficients",
             row_split = sub_grp, column_split = rs3,
             cluster_columns = F, 
             show_row_names = F,
             column_names_gp = gpar(fontsize = 5), 
             row_title_gp = gpar(fontsize = 9), 
             column_title_gp = gpar(fontsize = 8), 
             row_title_rot = 0,
             column_labels = col_labels[colnames(t(co_scale))],
             left_annotation = hleft, 
             cluster_rows = F, 
             # bottom_annotation = hbottom,
)
draw(hp)
dev.off()

