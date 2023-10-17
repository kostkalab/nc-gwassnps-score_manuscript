a2 = read_csv("./sup_data/sup_data_disease-snvs.csv.gz", col_types = cols(SNV_ID = col_character()))
a2$phenotype = tolower(gsub(" ","_",a2$phenotype))
## d1 = a2[a2$phenotype == "anorexia_nervosa",]
## d2 = a2[a2$phenotype == "autism_spectrum_disorder",]
## 
## intersect(d2$SNV_ID, d1$SNV_ID)
## length(intersect(d2$SNV_ID, d1$SNV_ID))
## duplicated(c(d2$SNV_ID, d1$SNV_ID))
## dim(a2)



library(ggplot2)
library(ggrepel)

# get coefficients and epigenome annotation
#------------------------------------------

coeff              <- readRDS("./data/analysis/coefficients.RDS")
top_mean           <- sapply(coeff, function(x) {x[,1]})
rownames(top_mean) <- rownames(coeff[[1]])
coeffi             <-  top_mean[-c(1,2),] # remove rowmeans and intercept, only concentrate on tissues


#- weights are bootstrap variances
wts           <- sapply(coeff, function(x) {x[,2]})
rownames(wts) <- rownames(coeff[[1]])
wts           <- wts[-c(1,2),] # remove rowmeans and intercept, only concentrate on tissues

library(erma) # https://bioconductor.org/packages/devel/bioc/vignettes/erma/inst/doc/erma.html
meta          <- mapmeta()
met           <- as.data.frame(meta[,c("Epigenome.ID..EID.", "GROUP","Standardized.Epigenome.name", "ANATOMY", "TYPE")])
rownames(met) <- met[,1]

# read in the file of GC correlation and model similarity

smat = readRDS("./data/analysis/coeffi_weighted_correlations.RDS")
sdat           <- data.frame(t(combn(rownames(smat),2)), stringsAsFactors = F) ### don't know why it doesn't work for me. have to change to charactor to work
sdat$sim       <- smat[cbind(sdat[,1],sdat[,2])] ### Nice way to do it!!
sdat           <- dplyr::tibble(sdat)
colnames(sdat) <- c("term_1","term_2","sim")
head(sdat)

gdat <- readRDS("./data/analysis/efo_gc.RDS")

strng_tls <- function(strng){
  #============================
  gsub("[[:punct:][:blank:]]+", " ", tolower(strng)) %>% stringr::str_replace_all("[[:space:]]+","_")
}

gdat <- gdat %>% mutate(term_name_id_1 = term_name_id_1 %>% tolower %>%  stringr::str_replace_all("[[:space:]]+","_"))
gdat <- gdat %>% mutate(term_name_id_2 = term_name_id_2 %>% tolower %>%  stringr::str_replace_all("[[:space:]]+","_"))
### sort alphabetically
gdat2 <- gdat %>% rowwise() %>% mutate(id_1 = min(term_name_id_1,term_name_id_2), id_2 = max(term_name_id_1,term_name_id_2))
sdat2 <- sdat %>% rowwise() %>% mutate(id_1 = min(term_1,term_2), id_2 = max(term_1,term_2))
smat_both <- inner_join(gdat2[,c(1,2,3,6,7)],sdat2[,3:5],by=c("id_1", "id_2"))

exmp = list(c("crohn's_disease", "ulcerative_colitis"), 
            c("anorexia_nervosa", "autism_spectrum_disorder"),
            c("obsessive-compulsive_disorder", "celiac_disease"))

# find the examples in smat_both
find_exmp = function(x){
  x = sort(x)
  which((smat_both$id_1 == x[1] & smat_both$id_2 == x[2])|(smat_both$id_2 == x[1] & smat_both$id_2 == x[2]))
}
id = sapply(exmp, find_exmp)





library(tidyverse)
plot_coeffi = function(X, wts, x1, x2, residual = T, n = 5, n0 = 3, pos, addin = NULL){
  md = quantile(wts,.25)/4
  twts =  1/sqrt(wts[     ,x1]*3/4+md)*1/sqrt(wts[     ,x2]*3/4+md)
  tmp <- lm(X[     ,x1] ~ X[     ,x2], weights =twts)
  # tmp <- lm(X[     ,x1] ~ X[     ,x2])
  
  df = data.frame(a1 = X[     ,x1], a2 = X[     ,x2], wt = twts)
  m = match(rownames(df), met$Epigenome.ID..EID.)
  df$name = paste(rownames(df), met$ANATOMY[m], sep = "-")
  
  # find the subset df that can represent something
  top_res = names(sort(abs(weighted.residuals(tmp)), decreasing = T)[1:5])
  subset_df = df[top_res,]
  summary(tmp)$r.squared
  
  # find the subset of df that contains the five lowest and five highest values
  
  subset_df = df[unique(c(order(df$a1)[c(1:n, (nrow(df)-n+1):nrow(df))], order(df$a2)[c(1:n, (nrow(df)-n+1):nrow(df))])),]
  if(nrow(subset_df) > 4*n0){
    t1 = order(df$a1)[1:n]
    t2 = order(df$a2)[1:n]
    l1 = order(df$a1)[127:(127-n+1)]
    l2 = order(df$a2)[127:(127-n+1)]
    
    # anything that appear twice will be picked out
    dup = c(t1,t2,l1,l2)[duplicated(c(t1,t2,l1,l2))]
    all = list(t1,t2,l1,l2)
    rest = lapply(all, function(x){
      names(x) = 1:n
      x1 = x[!x %in% dup]
      return(x1)
    })
    
    top = c(rest[[1]], rest[[2]])
    top = top[order(names(top))]
    lgh = n0*2 -floor((n*2-length(top))/2)
    top = top[1:lgh]
    low = c(rest[[3]], rest[[4]])
    low = low[order(names(low))]
    lgh = n0*2 -ceiling((n*2-length(low))/2)
    low = low[1:lgh]
    
    all = c(dup, top, low)
    print(length(all))
    subset_df = df[all,]
  }
  
  # add the addin tissue
  subset_df = rbind(subset_df, df[addin,])
  # add line
  exam = c(x1, x2)
  id = find_exmp(exam)
  sim = smat_both[id,]$sim
  rg = smat_both[id,]$rg
  reg_line = paste0("y=",
                    format(coef(tmp)[["X[, x2]"]], digits = 4),
                    "x + ",
                    format(coef(tmp)[["(Intercept)"]], digits = 3))
  textid = paste0("model_sim=",round(sim,2),"\n","genetic_corr=",round(rg,2),
                  "\n",reg_line)
  
  # ggplot(df,aes(x=a1, y=a2)) +
  #   geom_point() +
  #   geom_smooth(method='lm')
  
  p = ggplot(df, aes(x = a1, y = a2)) + 
    geom_point(aes(col=wt)) + 
    scale_colour_gradient(low = "grey89", high = "black", name = "")+ 
    geom_abline(slope = coef(tmp)[["X[, x2]"]], 
                intercept = coef(tmp)[["(Intercept)"]]) + 
    geom_text_repel(data = subset_df, 
                    aes(a1, a2, label=name),
                    size = 3.5)+
    labs(x = x1, y = x2) + 
    theme_minimal() + 
    theme(# The new stuff
      strip.text = element_text(size = 5), 
      axis.text = element_text( size = 10 ),
      axis.text.x = element_text( size = 10),
      axis.title = element_text( size = 14 ),
      legend.title=element_text(size=14), 
      legend.text=element_text(size=13)) 
  
  if (pos == "rightbottom"){
    p = p + annotate("text", x=Inf, y=-Inf, label= textid, hjust = 1, vjust = 0, size = 5)
  } else if(pos == "righttop"){
    p = p + annotate("text", x=Inf, y=Inf, label= textid, hjust = 1, vjust = 1, size = 5)
  } else if (pos == "lefttop"){
    p = p + annotate("text", x=-Inf, y=Inf, label= textid, hjust = 0, vjust = 1, size = 5)
  } else if (pos == "leftbottom"){
    p = p + annotate("text", x=-Inf, y=-Inf, label= textid, hjust = 0, vjust = 0, size = 5)
  }
  return(p)
}




pdf(file="./plot/model_sim_scatterplot.pdf", width = 6, height = 5)
p5 = plot_coeffi(X = coeffi, wts = wts, 
                 x1 = "crohn's_disease",
                 x2 = "ulcerative_colitis", 
                 pos = "lefttop")
p5

p2 = plot_coeffi(X = coeffi, wts = wts, 
                 x1 = "anorexia_nervosa",
                 x2 = "autism_spectrum_disorder", 
                 pos = "lefttop")

p2

p1 = plot_coeffi(X = coeffi, wts = wts, 
                 x2 = "celiac_disease",
                 x1 = "obsessive-compulsive_disorder", 
                 pos = "righttop", addin = "E116", n=4, n0=3)
p1




dev.off()

