mkdir data
cd data
wget ftp://ftp.sanger.ac.uk/pub/resources/software/gwava/v1.0/VEP_plugin/gwava_scores.bed.gz
wget ftp://ftp.sanger.ac.uk/pub/resources/software/gwava/v1.0/VEP_plugin/gwava_scores.bed.gz.tbi
# this score was used in DIVAN paper
# should be on hg19: see ftp://ftp.sanger.ac.uk/pub/resources/software/gwava/v1.0/README
# The results include the prediction scores from 3 different versions of the classifier, which are all in the range 0-1 with higher scores indicating variants predicted as more likely to be functional, and the underlying annotations used to compute these scores. For more details please refer to the GWAVA paper.






