library(ggplot2)

# GWAS catalog overall:
a2 <- readRDS("./data/GWAS_catalog_noncoding_hg19_uniquetrait.RDS") # read in the one that stores all gwas variants

catalog_ct100               <- readRDS("./data/GWAS_catalog_filtered_cutoff100_before_propagation.RDS") # read in the ones with out propagation >=100 snvs
catalog_all_snpsnap         <- readRDS("./data/catalog_cutoff100_snpsnap.RDS") # read all snpsnap matched snps

#### work on the object of the one before propagation. 
# find the snps that can be matched by snpsnap
catalog_1 = catalog_all_snpsnap[catalog_all_snpsnap$type == 1] # get only disease snps that can be matched by snpsnap
catalog_ct100_matched = lapply(catalog_ct100, function(x){
  x1 = x[x %in% catalog_1]
})
# only keep the diseases that have more than 100 snps
catalog_ct100_matched = catalog_ct100_matched[sapply(catalog_ct100_matched, length) >= 100]
names(catalog_ct100_matched) = NULL

catalog_ct100_matched = do.call("c", catalog_ct100_matched) # this answers how many SNPs are there without aggregation


#### Two objects, one before propagation and one after propagation
d      = catalog_ct100_matched
d_ag   = catalog_all_snpsnap[catalog_all_snpsnap$type == 1]

length(d) # number of snps non aggre ## 26080
length(d_ag) # number of snps aggre ## 77028

length(unique(d) ) # number of distinctive snps 20656
length(unique(d_ag)) # number of distictive snps 25516


## see the variants with the context
# i double checked on Jan 2022 that one snp annotated to two diseases will always has the same context such as intron or intergenic
table(d_ag %in% a2) # all true, say that all snps are in the GWAS catalog
# get the annotation
a3 = unique(a2[a2 %in% d_ag]) # all unique noncoding variants used in our study. 

context = a3$CONTEXT

context[context == "intergenic_variant x intron_variant"] = "intron_variant" # why? see line 55
context = gsub("acceptor_|region_|donor_", "", context) # 
context[context == "TF_binding_site_variant"] = "regulatory_variant"


df_context = data.frame(table(context))
names(df_context) = c("context", "number of variants")


df_context$percentage = df_context$`number of variants` / sum(df_context$`number of variants`)
library(scales)
df_context$percentage = percent(df_context$percentage, accuracy = 0.1)

df_context$context_nm = c("3'UTR", "5'UTR", "intergenic", "intron", "noncoding transcript", "regulatory region", "splice region")
# one intron X intergenic is actually intron


library(RColorBrewer)
library(cowplot)

df_context$percentage = df_context$`number of variants`/sum(df_context$`number of variants`)
df_context$context_nm = factor(df_context$context_nm, 
                               levels = c( "splice region", "5'UTR", "3'UTR", "noncoding transcript", 
                                           "regulatory region", "intergenic", "intron"))
df_context$label_y = c(0.9, 0.95, 0.75, 0.25, 0.85, 0.9,1.0)
df_context$label_x = c(2, 2, 1.1, 1.2, 2, 1,2)

df_context$x = 1.4
df_context$xend = 1.7
df_context$y = c(0.95,0.98,0.992,1, NA,NA,NA)
df_context$yend = c(0.85,0.9,0.95,1.0, NA,NA,NA)

p = ggplot(data=df_context, aes(x = NA, y=percentage, fill=context_nm)) +
  geom_bar(stat="identity") + 
  scale_fill_brewer(palette="Dark2", name = "Number of  \ndiseases associated") + 
  theme_classic() + 
  geom_text(aes(x = label_x-0.3, y = label_y, label=context_nm, hjust = 0), size=4) + 
  theme(legend.position = "none",
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) + 
  coord_cartesian(xlim = c(1,2)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  labs(title="", 
       x="GWAS catalog noncoding variants used in the study", y = "percentage") 



cowplot::save_plot(filename="./plot/stack_functional_class.pdf",p, base_height = 5, base_width = 4.5)

## write.csv(df_context,'~/project/variant-scores-2/data/GWAS_catalog_project_3.0/genetic/table1.csv', quote = F, row.names = F)

# before aggregation
d$snpID = paste(seqnames(d), start(d), sep = ":")
df1 = data.frame(table(table(d$snpID)))
names(df1) = c("number of diseases associated", "number of snps")

df2 = data.frame(table(table(d_ag$snpID)))
names(df2) = c("number of diseases associated", "number of snps")


gp1 = list(1,2:3,4:5,6:12)
gp2 = list(1,2:3,4:5,6:21)

get_sum = function(df, gp){
  dff = lapply(gp, function(x){
    y = df[df$`number of diseases associated` %in% x,]
    if (length(x) == 1){
      return(y)
    }else{
      y1 = data.frame(a1 = paste(min(x), max(x), sep = "-"), 
                      a2 = sum(y$`number of snps`))
      names(y1) = c("number of diseases associated", "number of snps")
      return(y1)
    }
  })
  dff = do.call("rbind", dff)
  dff$percentage = dff$`number of snps`/sum(dff$`number of snps`)
  return(dff)
}

dff1 = get_sum(df1, gp2)
dff2 = get_sum(df2, gp2)

## pie(dff1$`number of snps`, labels = dff1$`number of diseases associated`, 
##     main="Number of diseases associated with a SNV \n before propagation, (nSNVs = 20,656)")
## pie(dff2$`number of snps`, labels = dff2$`number of diseases associated`, 
##     main="Number of diseases associated with a SNV \n after propagation, (nSNVs = 25,516)")

# try barplot
## combine the dataframe
dff1$group = "before_propagation"
dff2$group = "after_propagation"

dff = rbind(dff1, dff2)
dff$group = factor(dff$group, levels = c("before_propagation", "after_propagation"))
dff$`number of diseases associated` = as.character(dff$`number of diseases associated`)
dff[dff$`number of diseases associated` == "6-21",]$`number of diseases associated` = ">5"
dff$`number of diseases associated` = factor(dff$`number of diseases associated`, 
                                             levels = c(">5", "4-5", "2-3","1"))

# writing
# less than 3 phenotypes. 
p2 = ggplot(data=dff, aes(x=group, y=`percentage`, fill=`number of diseases associated`)) +
  geom_bar(stat="identity") + 
  scale_fill_brewer(palette="YlGnBu", name = "Number of  \ndiseases associated") + 
  theme_classic() +
  labs(title="", 
       x="", y = "percentage") 
cowplot::save_plot(filename="./plot/stack_disease_numbers.pdf",p2, base_height = 5, base_width = 5)


# ggplot(data=dff, aes(x=`number of diseases associated`, y=`percentage`, fill=group)) +
#   geom_bar(stat="identity") + 
#   scale_fill_brewer(palette="Paired",
#                     name = "group", labels = c("before propagation, n=20,656",
#                                                "after propagation, n=26,080")) + 
#   theme_minimal() +
#   labs(title="", 
#        x="number of diseases associated", y = "percentage")+
#   theme(legend.position = c(0.75, 0.6))
# 

