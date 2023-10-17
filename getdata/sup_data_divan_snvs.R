library(readr)

combine_41_snpsnap <- readRDS("~/project/variant-scores-2/data/DIVAN/data/combine_41_snpsnap.RDS")
combine_41_snpsnap <- readRDS("./data/divan_41_snpsnap.RDS")

# sort by trait name
nms = as.character(sapply(combine_41_snpsnap, function(x) x$trait[1]))
combine_41_snpsnap = combine_41_snpsnap[order(nms)]

names(combine_41_snpsnap) = NULL
gr_divan = do.call("c",combine_41_snpsnap)
tmp = gr_divan[gr_divan$type == 1]
t = tapply(as.factor(tmp$divan_new), tmp$trait, table)
t = do.call("rbind", t)
t = t[t[,"0"] >= 50 & t[,"1"] > 20, ]

gr_divan = gr_divan[gr_divan$trait %in% rownames(t)]
gr_divan$divan_new[gr_divan$divan_new == 0] = "train"
gr_divan$divan_new[gr_divan$divan_new == 1] = "test"


df = data.frame(SNV_ID = gr_divan$snpID,
                rsID = gr_divan$rsID,
                phenotype = gr_divan$trait,
                hg19_chromosome = seqnames(gr_divan),
                hg19_location = start(gr_divan), 
                type = gr_divan$type,
                input_snp = gr_divan$input_snp,
                context = gr_divan$divan_new)

colnames(df)[7] = "disease_SNV_ID"
colnames(df)[6] = "disease_associated"
df$disease_associated = as.logical(df$disease_associated) 
write_csv(df ,'./sup_data/sup_data_divan-snvs.csv.gz')

