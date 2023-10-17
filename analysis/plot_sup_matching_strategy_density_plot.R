library(ggplot2)

g1kg = readRDS("./data/1KG_processed_eur.RDS")

tss_19 = readRDS("./data/tss_19.RDS")
ld0.8 = read.table(gzfile("./snpsnap/ld0.8_collection.tab.gz"), sep = "\t", header = T)

tss        = readRDS("./data/catalog_control_tss.RDS")
snpsnaptss = readRDS("./data/catalog_control_snpsnap_tss.RDS")
random     = readRDS("./data/catalog_control_random.RDS")

# get the distance to nearest protein coding genes by using TSS 19 i generated before
dist      = distanceToNearest(snpsnaptss, tss_19)
snpsnaptss$dist = mcols(dist)[,1]

dist       = distanceToNearest(random, tss_19)
random$dist = mcols(dist)[,1]

# get the distance to nearest genes in snpsnap by using ld0.8 database
tss$dist_gene = ld0.8[match(tss$snpID, ld0.8$snpID),"dist_nearest_gene"]
snpsnaptss$dist_gene = ld0.8[match(snpsnaptss$snpID, ld0.8$snpID),"dist_nearest_gene"]
random$dist_gene = ld0.8[match(random$snpID, ld0.8$snpID),"dist_nearest_gene"]

dis   = tss[tss$type == 1]
c_tss = tss[tss$type == 0]
c_snpsnaptss = snpsnaptss[snpsnaptss$type == 0]
c_random     = random[random$type == 0]

## now plot the distance to nearest protein coding genes
dat1 = data.frame(cond = "disease", rating = dis$dist)
dat2 = data.frame(cond = "control: tss matched", rating = c_tss$dist)
dat3 = data.frame(cond = "control: snpsnap_tss matched", rating = c_snpsnaptss$dist)
dat4 = data.frame(cond = "control: random matched", rating = c_random$dist)

dat = rbind(dat1, dat2, dat3, dat4)
dat$cond = factor(dat$cond, levels = c("disease", "control: tss matched", "control: snpsnap_tss matched",
                                       "control: random matched"))
# Density plots with means
p1 = ggplot(dat, aes(x=rating, colour=cond)) +
  geom_density() +
  scale_x_continuous(trans='log10') + 
  labs(x = "distance to the nearest TSS(bp), protein-coding genes", 
       y = "density of SNVs",
       colour = "matching strategy")


## now plot the distance to nearest all genes in snpsnap
dat1 = data.frame(cond = "disease", rating = dis$dist_gene)
dat2 = data.frame(cond = "control: tss matched", rating = c_tss$dist_gene)
dat3 = data.frame(cond = "control: snpsnap_tss matched", rating = c_snpsnaptss$dist_gene)
dat4 = data.frame(cond = "control: random matched", rating = c_random$dist_gene)

dat = rbind(dat1, dat2, dat3, dat4)
dat$cond = factor(dat$cond, levels = c("disease", "control: tss matched", "control: snpsnap_tss matched",
                                       "control: random matched"))
# Density plots with means
p2 = ggplot(dat, aes(x=rating, colour=cond)) +
  geom_density() +
  scale_x_continuous(trans='log10') + 
  labs(x = "distance to the nearest TSS(bp), all genes", 
       y = "density of SNVs",
       colour = "matching strategy")
# there are 18187 rows removed because not all tss and random control variants are in the ld0.8 database. so na is generated. 


pdf("./plot/matching_strategy_density_plot.pdf", width = 7, height = 5)
p1
p2
dev.off()
