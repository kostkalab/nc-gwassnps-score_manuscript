library(readr)
library(dplyr)
library(readxl)

a15 = read_csv("./sup_data/sup_data_pairwise-divan-individual.csv.gz")

aa = dplyr::filter(a15, score1 == "DIVAN" & score2 == "GenoCanyon" & measure == "pr")
aa2 = dplyr::filter(aa, higher_median == "GenoCanyon")

aa2 = aa %>% select(disease, higher_median) %>% mutate(GenoCanyon_better_our = (higher_median == "GenoCanyon"))

library("readxl")
divan_data <- read_excel("./external_data/13059_2016_1112_MOESM1_ESM.xlsx",
                      sheet = "TableS16")
colnames(divan_data)[10] = "diff"

# are they using the same name? 
# i think so, because efo term, lupus is called system lupus erythematosus

divan_data$disease = tolower(gsub(" ", "_", divan_data$disease))

divan_data2 = divan_data %>% select(disease, DIVAN, GenoCanyon) %>% 
  mutate(GenoCanyon_better_chen = (GenoCanyon > DIVAN))

# combine them together
df_all = full_join(aa2, divan_data2, by = "disease") %>% select(disease, GenoCanyon_better_our, GenoCanyon_better_chen)

df_all = inner_join(aa2, divan_data2, by = "disease") %>% 
  select(disease, GenoCanyon_better_our, GenoCanyon_better_chen) %>% 
  mutate(GenoCanyon_better_both = GenoCanyon_better_our & GenoCanyon_better_chen)

df_all = inner_join(aa2, divan_data2, by = "disease") %>% 
  select(disease, GenoCanyon_better_our, GenoCanyon_better_chen) 

write_csv(df_all, file = "./sup_data/sup_data_perf-divan-our-vs-chen.csv.gz")
