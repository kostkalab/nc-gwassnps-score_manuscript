smat = readRDS("./data/analysis/coeffi_weighted_correlations.RDS")
sdat           <- data.frame(t(combn(rownames(smat),2)), stringsAsFactors = F) ### don't know why it doesn't work for me. have to change to charactor to work
sdat$sim       <- smat[cbind(sdat[,1],sdat[,2])] ### Nice way to do it!!
sdat           <- dplyr::tibble(sdat)
colnames(sdat) <- c("term_1","term_2","sim")
head(sdat)

gdat <- readRDS("./data/analysis/efo_gc.RDS")
head(gdat)

library(readr)
a4 = read_csv("./sup_data/tab_model-vs-gc.csv")


library(dplyr)

strng_tls <- function(strng){
  #============================
  gsub("[[:punct:][:blank:]]+", " ", tolower(strng)) %>% stringr::str_replace_all("[[:space:]]+","_")
}

gdat <- gdat %>% mutate(term_name_id_1 = term_name_id_1 %>% tolower %>%  stringr::str_replace_all("[[:space:]]+","_"))
gdat <- gdat %>% mutate(term_name_id_2 = term_name_id_2 %>% tolower %>%  stringr::str_replace_all("[[:space:]]+","_"))

gdat %>% head()

# smat_both <- inner_join(gdat,sdat,by=c("term_name_id_1"  ="term_1", "term_name_id_2"="term_2")) ### there's something wrong with this one
### i add these
dim(gdat)
dim(sdat)
### end

### sort alphabetically
gdat2 <- gdat %>% rowwise() %>% mutate(id_1 = min(term_name_id_1,term_name_id_2), id_2 = max(term_name_id_1,term_name_id_2))
sdat2 <- sdat %>% rowwise() %>% mutate(id_1 = min(term_1,term_2), id_2 = max(term_1,term_2))

smat_both <- inner_join(gdat2[,c(1,2,3,6,7)],sdat2[,3:5],by=c("id_1", "id_2"))
exmp = list(c("crohn's_disease", "ulcerative_colitis"), 
            c("anorexia_nervosa", "autism_spectrum_disorder"),
            c("celiac_disease", "obsessive-compulsive_disorder"))

# find the examples in smat_both
find_exmp = function(x){
  which((smat_both$id_1 == x[1] & smat_both$id_2 == x[2])|(smat_both$id_2 == x[1] & smat_both$id_2 == x[2]))
}

tab_examp = apply(a4, 1, function(x) x[1:2])
id2 = apply(tab_examp, 2,find_exmp)

id = sapply(exmp, find_exmp)
col = rep(rgb(0,0,0,1/2), nrow(smat_both))
col[id2] = "#A13D2D"
col[id2] = rgb(1,0,0,1/2)
col[id] = rgb(1,0,0)
cex = rep(.75,nrow(smat_both))
cex[id] = 1.25

pdf(file="./plot/mol_sim_gc_comparison2.pdf", width = 6, height = 6)
plot(smat_both$rg, smat_both$sim, xlim=c(-1/2,1),ylim=c(-1/2,1),
     cex=cex,pch=19,col=col, 
     xlab="genetic correlation",ylab="model correlation",
     cex.axis=1,
     cex.lab=1.5,
     mgp=c(2.5,1,0));
points(smat_both$rg[id], smat_both$sim[id], 
       cex=.75,pch=19,col=rgb(0,0,0,1/2))
grid()
abline(a=0,b=1,col=2)
abline(h=0.27,lt=2)
abline(v=0.19,lt=2)
text(smat_both$rg[id]-0.07, smat_both$sim[id], labels = c("B", "C", "D"), font=2, cex = 1.5, col = "blue")
# add rectangle 
# rect(0.25, -0.5, 1, 0.25,
#      col= rgb(red = 0.95, green = 0.5, blue = 0, alpha = 0.2)) # right bottom
# 
# rect(0.25, 0.25, 1, 1,
#      col= rgb(red = 0.1, green = 0.1, blue = 0.9, alpha = 0.2)) # right up blue
# 
# rect(-0.5, 0.25, 0.25, 1.0,
#      col= rgb(red = 0.1, green = 0.9, blue = 0.1, alpha = 0.2)) # right left green

dev.off()


cor.test(smat_both$rg, smat_both$sim) #- r=0.32 p-value = 2.416e-15
cor.test(smat_both$rg, smat_both$sim, method="spearman") #- 0.14 p-value = 0.0006388
# inteprete: https://bookdown.org/ndphillips/YaRrr/correlation-cor-test.html
# https://rcompanion.org/rcompanion/e_01.html

model = lm(rg ~ sim,
           data = smat_both)

summary(model) 

