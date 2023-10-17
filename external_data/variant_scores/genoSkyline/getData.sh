
#- According to the http://genocanyon.med.yale.edu/GenoSkyline these genosSkyLinePlus scores are 
#
#  - hg19
#
#- and from looking at the files we infer they are
#
#  - zero-based, half-open [,) 
# 



# 1) get data

wget -nc http://genocanyon.med.yale.edu/GenoSkylineFiles/GenoSkylinePlus/GenoSkylinePlus_bed.tar.gz

# 2) unzip 

gunzip GenoSkylinePlus_bed.tar.gz

# 3) extract

tar xvf GenoSkylinePlus_bed.tar

# 4) bgzip (for tabbix later)

cd BedGraph
find . -name "*.bedGraph" -print0 | xargs -0 -I {}  sh -c "grep -v -E '^#|^$' {} | sort -k1,1 -k2,2n -k3,3n | bgzip > {}.gz"
find . -name "*.bedGraph" -print0 | xargs -0 -I {}  rm {}
 
#- 5) tabbix ; assuming zero-base half-open because the files start with zero and repeat the end.

find . -name "*.bedGraph.gz" -print0 | xargs -0 -I {} tabix -0 -b 2 -e 3 -s 1 {}

#- clean up

cd ..
rm GenoSkylinePlus_bed.tar



