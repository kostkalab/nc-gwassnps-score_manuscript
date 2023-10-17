# README
------------------------------------------

This repository stores the code to get the data, plot and results for the paper: Disease-specific analysis improves prioritization of non-coding genetic variants.

## Objectives ##
 
In the paper *Disease-specific analysis improves prioritization of non-coding genetic variants*, we presents a strategy that uses tissue-specific scores
to generate disease-specific variant scores that can be used to prioritize variants for a specific disease you are interested in. 

In this repository, we share the code and data that can do the following: 

* Generate disease-associated and control SNVs for 111 diseases.
* Run the regularized logistic regression model to learn and predict disease-associated variants for those diseases.
* Analyze disease-specific model.


## Environment ##

You can get the code from bitbucket by running the following command:

```
$ git clone https://QianqianLiang@bitbucket.org/dkostka/variant-scores_manuscript.git
$ cd scds_manuscript

```

All the linux command and R code are under this current directory. 

If you want to use our processed data, you can find the data in the folder `./data/` on xxx website.

If you want to use the original data, such as the GWAS catalog, 1000 genome project, different variant scores, etc., you can find the code to get the data under `./external_data/`. 


## Use ##

### Overview ###

This is an overview of the structure of the repository. 

The following folders store the code: 

* `./getdata`   -> generate the disease associated and control variants for 111 diseases.
* `./model`     -> run the regularized logistic regression models on variants for 111 diseases. 
* `./analysis`  -> analyze the beta coefficients of the model. 

The following folder stores the data and software:

* `./external_data`  -> external data used
* `./snpsnap`    -> snpsnap software

The following folders store our processed dataset and plots:

* `./data`       -> processed data we generated in every steps. 
* `./sup_data`   -> tables and supplemental tables that are published in the paper
* `./plot`       -> original plots published in the paper.

All the data is from ./data and ./sup_data folder is in zenodo. This is an temperary url: <10.5281/zenodo.8191149>


### Step 1: Get disease and control SNVs ###

#### External data annotation ####



#### Diesase SNVs ####

We applied the following steps on **the GWAS Catalog** data to get the disease-associated SNVs in my study: 

* Filter for SNVs (deletions, duplications, etc. was excluded). - stored as `./data/GWAS_catalog_hg19.RDS` 
* Filter for noncoding SNVs. - stored as `./data/GWAS_catalog_noncoding_hg19.RDS`
* Annotate one SNV to a single efo term. (some SNVs are annotated to a trait with more than one efo terms.) - stored as `./data/GWAS_catalog_noncoding_hg19_uniquetrait.RDS`
* Filter only for EFO traits (trait ID with EFO:xxxxxxxxx).
* Propagate traits to include descendants. - stored as `./data/GWAS_catalog_noncoding_hg19_uniquetrait_propagate.RDS`
* Filter only for a 'disease' term. 
* keep diseases with >100 SNVs - stored as `./data/GWAS_catalog_filtered_cutoff100.RDS`
* keep diseases with >100 SNVs before propagation - stored as `./data/GWAS_catalog_filtered_cutoff100_before_propagation.RDS`

The code are stored as:

* `./getdata/get_noncoding_hg19.R`
* `./getdata/get_GWAS_catalog_hg19.R` -> this generate `GWAS_catalog_hg19.RDS` and `GWAS_catalog_noncoding_hg19.RDS`
* `./getdata/get_GWAS_catalog_hg19_uniquetrait.R`
* `./getdata/get_GWAS_catalog_noncoding_hg19_uniquetrait_propagate.R`
* `./getdata/get_GWAS_catalog_filtered_cutoff100.R`




#### Control SNVs ####


For each SNV, we generate 10 control variants using four different matching strategies, sepearately: 
 
 
1. SNPsnap_TSS
2. SNPsnap
3. TSS
4. Random

For 1 and 2, we used a re-implementation of [SNPsnap](https://data.broadinstitute.org/mpg/snpsnap/match_snps.html). 
The code to run SNPsnap is under the folder `./snpsnap/` 

For 3 and 4, we got the matched control noncoding SNVs from **1000 Genome Project**, matched with protein-coding TSS or random. 


##### SNPSNAP #####

###### 1. Run SNPsnap software ###### 

To run SNPsnap, download the original LD annotation file from the SNPsnap website [link][https://data.broadinstitute.org/mpg/snpsnap/database_download.html]. 
Here I used the LD0.8 European as threshold for distance cut off(`ld0.8_collection.tab.gz`), but you can use other threholds as well. 

Input file format is `chr:pos` such as `1:124987986`, while each SNV takes a line. The first line of the file is `lead_snps`. 

Use `run_match_snps.sh` to specify the position of your input file, output file, annotation file etc. 

If you want to change the threshold of the four criteria, i.e. distance to the nearest gene, gene density, LD buddies and minor allele frequency, 
change the line 294 of the python script `match_snps.py`: 
```
thresholds = set_thresholds(max_dev_maf=0.05, max_dev_gene_count=0.50, max_dev_dist_nearest_gene=0.20, max_dev_friends_ld=0.50)
```

Please see here [link][https://data.broadinstitute.org/mpg/snpsnap/documentation.html] for the meaning of those four criteria. 

To run SNPsnap, run the following command: 

```

conda create --name SNIPSNAP python=3
conda activate SNIPSNAP
pip install numpy
pip install pandas
cd ./snpsnap/v2_2020_04_27/
chmod u+x ./run_match_snps.sh
./run_match_snps.sh

```

In my study, for **SNPsnap_TSS**, i set the threshold as following: 
```
thresholds = set_thresholds(max_dev_maf=5, max_dev_gene_count=10000, max_dev_dist_nearest_gene=0.20, max_dev_friends_ld=10000)

```

For **SNPsnap**, I set the threshold as following: 
```
thresholds = set_thresholds(max_dev_maf=0.05, max_dev_gene_count=0.50, max_dev_dist_nearest_gene=0.20, max_dev_friends_ld=0.50)

```

For each SNP, I set the program to match with 100 SNVs (if less than 100 SNVs match the requirement, there will be duplicates). 
I will later select 10 matched control SNVs from them. 


###### 2. Annotate input and matched control SNVs ###### 

Next, annotate disease and matched control SNVs using the file `./snpsnap/ld0.8_collection.tab.gz`. 
The annotated files are under: `./data/snpsnap_annotated`

The code to annoate SNVs is 

* `./getdata/fun_snpsnap_annotation.R`
* `./getdata/get_snpsnap_annotated.R`

###### 3. Randomly Select 10 matched control SNVs ###### 

After this, I matched each input SNV with 10 snpsnap matched SNVs. 
The data is stored as: 

* `./data/catalog_control_snpsnap.RDS`
* `./data/catalog_control_snpsnaptss.RDS`

The code to get 10 matched control is: 

* `./getdata/get_catalog_control_snpsnap.R`


##### 1000 Genome #####

For 3 and 4, I get matched control variants from SNVs in 1000 Genome project. 
I applied the following filters for the 1000 Genome SNVs: 

* is a SNV
* minor allele frequency in european population >0.01 -> stored as `./data/1KG_processed_eur.RDS`
* is a noncoding SNV -> stored as `./data/1KG_processed_eur_noncoding.RDS`

In the files above, each SNV is annotated with the distance to the nearest transcription start site (TSS) nearest protein-coding gene. 

You can find transcription start site data under: 
* `./data/TSS_19.RDS` 
* `./data/TSS_38.RDS`

Here, I used gencode to get transcripts and used ensemble to get protein-coding transcripts. The code stored as `./getdata/get_TSS19_and_TSS38.R`
Since the ensemble dataset keeps updating, if you want to replicate the same result of the distance of nearest TSS, use the above TSS file we provided directly. 

For **3.tss** matching, we match each SNV to 10 SNVs that has similar distribution to nearest TSS. 

For **4.random** matching, we match each SNV to 10 10 randomly selected SNVs from filtered 1000 Genome SNVs. 

The code to generate them are under: 

* `./getdata/get_catalog_control_tss.R`
* `./getdata/get_catalog_control_random.R `

##### disease and control variants for 111 diseases #####

Now we have 10 matched SNVs for each disease-associated SNV in four strategies. 
Then we will retrieve all disease and control SNVs for 111 diseases. 

* `./data/catalog_cutoff100_snpsnap.RDS`
* `./data/catalog_cutoff100_snpsnaptss.RDS`
* `./data/catalog_cutoff100_tss.RDS`
* `./data/catalog_cutoff100_random.RDS`

The code to get them `./getdata/get_catalog_cutoff100_four_controls_bydisease.R`


#### LD block annotation ####

For each disease-associated SNV, we annotate the LD block id the SNV resides in, p.value reported in the GWAS catalog, 
and whether it is a 'representative SNV' (SNV with the smallest pvalue in the LD block). 

The control SNVs will be annotated the same as the matched disease-associated SNV. 

We generate the data: 

* `./data/catalog_cutoff100_snpsnap_ld_marked.RDS`

The code that generate the data: 

* `/getdata/fun_ld_marked`
* `/getdata/get_catalog_cutoff100_snpsnap_ld_markded.R`


#### Chromosome heldout ####

In the chromosomal heldout setting, training and test SNVs are not in the same chromosome. 
Here, we annotate each SNV whether it is in the training or test set.
For each disease, we select a group of chromosomes that has around 20% SNVs with a 1/10 positive and negative ratio (the same as in cross-validation setting).
We marked rest of the SNPs as the training set.

We generate the data: 

* `./data/catalog_cutoff100_snpsnap_chrom_heldout.RDS`

The code that generate the data: 

* `/getdata/fun_chrom_heldout_annotation`
* `/getdata/get_catalog_cutoff100_snpsnap_chromosome_heldout.R`

#### DIVAN dataset ####

The DIVAN dataset contains test and control variants for 41 traits used to compare our model with DIVAN publication. (I used 29 disease traits with more than 10 test and 50 train diseases SNPs for further analysis.)

To generate the DIVAN dataset, I used two external datasets and one mapping dataset: 

* ARB datasest from [Association Results Browser](https://www.ncbi.nlm.nih.gov/projects/gapplus/sgap_plus.htm). 
* phegen dataset from [Phenotype Genotype Integrater](https://www.ncbi.nlm.nih.gov/gap/phegeni)
* mapping dataset from [OxO](https://www.ebi.ac.uk/spot/oxo/)

I processed them to be used later by DIVAN, generating files: 

* `/data/ARB_processed.RDS` <- generated by `./getdata/get_ARB_processed.R`
* `./data/phegen_processed.RDS` <- generated by `./getdata/get_phegen_processed.R`
* `./data/mesh_efo_mapping_45.RDS` <- generated bu `./getdata/get_mesh_efo_mapping_45.R`

I used them to generate the disease-associated variants and matched them with snpsnap (see previous section for snpsnap matching), generating the follow data: 

* `./data/divan_41.RDS` <- generated by `./getdata/get_DIVAN_41.R`
* `./data/divan_41_snpsnap.RDS` <- generated by `./getdata/get_DIVAN_snpsnap_annotated.R` and `./getdata/get_DIVAN_41_snpsnap.R`

I then annotate them with organism-level, tissue-specific and DIVAN score, generating the following data: 

* `./data/anno_divan_DHS.RDS` 
* `./data/anno_divan_DIVAN.RDS`
* `./data/anno_divan_GenoCanyon.RDS` <- all three generated by `./getdata/anno_divan.R`


### Step 2: Get annotation for SNVs ###

#### Organism-level score annotation ####

We provide annotation for five organism-level variant scores:

* CADD
* eigen
* GWAVA
* GenoCanyon
* LINSIGHT

#### Tissue-specific score annotation ####

We provide annotation for three tissue-specific variant scores: 

* DHS (DNase hypersensitive sites)
* Genoskyline
* Fitcons2

All the precomputated organism and tissue-specific scores are under `./external_data/variant_scores` and `./external_data/avocado`

#### Disease-specific score annotation ####

We also get the precomputed disease-specific variant score:

* DIVAN

#### Code and data ####

Code that can annotate SNVs with organism-level and tissue-specific scores are at: 

* `./getdata/fun_read_variantscore_organism.R`
* `./getdata/fun_read_variantscore_tissue.R`
* `./getdata/fun_read_DIVAN.R`
* `./getdata/anno_catalog_cutoff100_orgamism.R`
* `./getdata/anno_catalog_cutoff100_tissue.R`

Annotated SNVs for disease-associated SNVs and snpsnap matched control SNVs are stored at: 

* `./data/anno_catalog_cutoff100_snpsnap_orgamism.RDS`
* `./data/anno_catalog_cutoff100_snpsnaptss_orgamism.RDS`
* `./data/anno_catalog_cutoff100_tss_orgamism.RDS`
* `./data/anno_catalog_cutoff100_random_orgamism.RDS`
* `./data/anno_catalog_cutoff100_snpsnap_DHS.RDS`
* `./data/anno_catalog_cutoff100_snpsnap_genoskyline.RDS`
* `./data/anno_catalog_cutoff100_snpsnap_Fitcons2.RDS`


### Step 3: Organism-level variant scores performance ###

#### Bootstrap samples #### 

To compare different variant scores, we generate 30 bootstrap samples for each disease to measure the performance of variant scores. 
Each bootstrap sample randomly contains 90% of disease and control variants. 

The following code generate bootstrap samples: 

* `./analysis/get_bootstrap_samples.R` -> generating file `./data/analysis/bootstrap_samples_30_0.9.RDS`


#### Organism-level variant score performance ####

The AUC of different organism-level scores on 111 diseases is under: 

* `./data/analysis/perf_org.RDS`
* `./data/analysis/bootstrap_org_auc.RDS` - organism-level score performance for boostrap samples

The code to get them is under: 

* `./analysis/fun_perf_org`
* `./analysis/perf_org`
* `./analysis/bootstrap_org_auc.R`



### Step 4: Build the model for tissue-specific variant scores: baseline and tissue-weighted (disease-specific) ###

We build a regularized logistic regression model (named *tissue-weighted*) to learn and predict disease-associated variants for a specific disease. 
Overall we have the following different settings: 

To train and test the model, we use cross-validation. We have two settings here:

1. Cross-validation. -`./model/catalog_cv.R`
2. Cross-validation with variants not in the same LD block -`./model/catalog_cv_ldfree.R`
7. Cross-validation, measure the performance in bootstrap samples/ -`./model/catalog_cv_bootstrap.R`

We also train and test the model using test SNVs not in the same chromosome of the training SNVs: 
3. Chromosome held-out. (variants for train and test are on different chromosome) -`./model/catalog_chr_holdout.R`
4. Under chromosome held-out setting, use randomly selected SNVs as control. 

We also did: 

5. Use all variants to train the model and get prediction scores. -`./model/catalog_predscore.R` 
6. Get the coefficients of the model.  -`./model/catalog_get_coefficient.R`


The data that are generated are: 

* `./data/analysis/perf_DHS_logreg.RDS` - cross-validation setting
* `./data/analysis/perf_Fit_logreg.RDS`
* `./data/analysis/perf_geno_logreg.RDS`
* `./data/analysis/perf_DHS_logreg_ldfree.RDS` - cross-validation, with variants not in the same LD block
* `./data/analysis/perf_Fit_logreg_ldfree.RDS`
* `./data/analysis/perf_geno_logreg_ldfree.RDS`
* `./data/analysis/bootstrap_DHS_SNPSNAP_logreg.RDS` - tissue-weighted performance in bootstrap samples
* `./data/analysis/bootstrap_Fit_SNPSNAP_logreg.RDS`
* `./data/analysis/bootstrap_geno_SNPSNAP_logreg.RDS`
* `./data/analysis/perf_DHS_chrom.RDS` - chromosome heldout setting. 
* `./data/analysis/perf_DHS_random.RDS` - randomly select test SNVs, under chromosome-heldout setting
* `./data/analysis/model_111diseases_DHS.RDS` - our model for 111 diseases
* `./data/analysis/model_111diseases_Fitcons2.RDS`
* `./data/analysis/model_111diseases_Genoskyline.RDS`
* `./sup_data/sup_data_prediction-scores-dhs-weighted.csv.gz` - prediction scores
* `./data/abalysis/coefficients.RDS` - coefficients of our model, including standard deviation.


### Step 5: Analysis ###

#### Pairwise comparison ####

To compare different variant scores, we test their performances on bootstrap variant samples. 
Each bootstrap contains 90% of all variants (90% for positive variants, and 90% for control variants).
We use this to test the performance for organism-level and tissue-specific scores. 

We provide pairwise comparison for the following senerios: 

* Org vs. org
* Tiss vs. Tiss
* Tiss vs. Org
* DIVAN vs. GenoCanyon vs Tissue-weighted-DHS (independent dataaset)


The code to generate pairwise comparison is under: 

* `./analysis/pairwise_divan.R`
* `./analysis/pairwise_org_aggre.R`
* `./analysis/pairwise_org_indiv.R`
* `./analysis/pairwise_tis_aggre.R`
* `./analysis/pairwise_tis_indiv.R`
* `./analysis/pairwise_tis_vs_org_aggre.R`
* `./analysis/pairwise_tis_vs_org_indiv.R`

The pairwise comparison results is under 

* `./supdata/sup_data_pairwise-*`

#### DIVAN performance ####

* To generate bootstrap sample for DIVAN performance: `./analysis/divan_get_bootstrap_samples.R` -> `./data/analysis/bootstrap_samples_divan29.RDS`
* To get the performance of DIVAN vs GenoCanyon vs DHS-weighted and pairwise comparison statistics: `./analysis/divan_get_auc.R` -> `./data/analysis/p_value_DIVAN_top.RDS`


#### Coefficient analysis ####

In this section, we analyzed the coefficients that are derived from the regularized logistic regression model in part 3. 
We get the weighted correlations for the coefficients, and then use them to make diseases into seven clusters. 
For each cluster, we get top 5 tissues that can distinguish this cluster from other clusters. 
In the last, we also the genetic correlations of those 111 diseases and compare them with the model similarity. 

All the code and data that we generated are as below: 

* To get the weighted correlations(model similarities): `./analysis/coeffi_weighted_correlations.R` -> `./data/analysis/coeffi_weighted_correlations.RDS`
* To get the seven clusters: `./analysis/coeffi_clusters.R` -> `./data/analysis/coeffi_clusters.RDS`
* To get a name for seven clusters: `./analysis/coeffi_cluster_name.R` -> `./supdata/sup_data_cluster_term-frequency.csv.gz`
* To get the top tissues: `./analysis/coeffi_top_tissues.R` -> `./data/analysis/coeffi_top5tissues.RDS`
* To get the genetic correlation: `./analysis/gc_genetic_correlation.R` -> `./data/analysis/efo_gc.RDS`


### Plot ###

All plots for the paper are under `./plot/`. The code that generates the plot are under `./analysis/` with a matching name start with `plot_`. 

### Table, supplemental table and  supplemental data ###

All the table, supplemental table and supplemental data are under the folder: `./sup_data/`

The code that generate those data are:

* `./getdata/sup_data_*` 
* `./analysis/sup_data_*`
* `./analysis/pairwise_divan.R`
* `./analysis/pairwise_org_aggre.R`
* `./analysis/pairwise_org_indiv.R`
* `./analysis/pairwise_tis_aggre.R`
* `./analysis/pairwise_tis_indiv.R`
* `./analysis/pairwise_tis_vs_org_aggre.R`
* `./analysis/pairwise_tis_vs_org_indiv.R`
* `./analysis/coeffi_cluster_name.R`
* `./model/catalog_predscore.R`

