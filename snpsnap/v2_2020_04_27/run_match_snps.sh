#!/bin/bash


input_snps_file="input_snps.txt"
n_matches=10
ldbud_r2="friends_ld05"
output_root="output"
analysis_name="example"

db_file="ld0.1_collection.tab.gz"
# echo "python match_snps.py $input_snps_file $n_matches $ldbud_r2 $output_root $analysis_name"

python match_snps.py $input_snps_file $n_matches $ldbud_r2 $db_file $analysis_name $output_root

