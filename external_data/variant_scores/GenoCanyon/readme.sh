mkdir data

cd data

wget http://genocanyon.med.yale.edu/GenoCanyonFiles/GenoCanyon_Chr{1.22}_hg19.tar.gz
wget http://genocanyon.med.yale.edu/GenoCanyonFiles/GenoCanyon_Chr{X,Y}_hg19.tar.gz

# from website: http://genocanyon.med.yale.edu/GenoCanyon_Downloads.html

# make tabix file
for file in $(ls ./*/*.tsv); do tr ' ' '\t' < $file | cat > ${file%%.txt}.tsv ; done
for file in $(ls ./*/*.tsv); do bgzip $file ; echo $file; done
for file in $(ls ./*/*.tsv.gz); do tabix -s 1 -b 2 -e 2 $file ; echo $file; done





