#!/bin/bash


## This script takes multiple VCF files as input, extracts SNV and INDELS variants, merges variants, and generates a unified mutation table as output.
## The output can be used as an input for analysis such a non-synonymous/synonymous mutation cancer driver analysis, e.g., dndscv (Martincorena et al., 2017).

#######################################

# INDEL ANALYSIS

# step1 
for file in *vcf.gz
do

echo ${file} >> vcf_files_list.txt

done

#step2

cat vcf_files_list.txt | while read line;
do
bcftools view -v snps ${line} -Oz -o ${line:0:4}_snp.vcf.gz;
bcftools view -v indels ${line} -Oz -o ${line:0:4}_indels.vcf.gz;
done

#step3

for file in *indels.vcf.gz
do

echo ${file} >> indels_list.txt
done

cat indels_list.txt | while read line;
do

sampleId=${line:0:4}

zcat $line | grep -v "##" | awk '{gsub (/#/,"")}1' | awk -v ID="${sampleId^}" 'NR==1{print "SampleID", $1, $2, $4, $5; next}{print ID, $1, $2, $4, $5}' | tail -n +2 > ${sampleId}.indels.txt
done

#concatenates all *_final_indel files and attach the dndscv headers 

for file in *indels.txt ; do cat $file >> combined.txt ; done
awk 'BEGIN {print "SampleID chr pos ref mut"}{print $0}' combined.txt  > indel_mutation_table.txt
rm combined.txt
rm indels_list.txt
rm *_indels.vcf.gz
rm *.indels.txt


# SNP ANALYSIS


#step4


for file in *snp.vcf.gz
do

echo ${file} >> snp_list.txt
done

cat snp_list.txt | while read line;
do

sampleId=${line:0:4}

zcat $line | grep -v "##" | awk '{gsub (/#/,"")}1' | awk -v ID="${sampleId^}" 'NR==1{print "SampleID", $1, $2, $4, $5; next}{print ID, $1, $2, $4, $5}' | tail -n +2 > ${sampleId}.snp.txt
done
# concatenates all *_final_indel files and attach the dndscv headers 

for file in *snp.txt ; do cat $file >> snp_mutation_table.txt ; done

rm snp_list.txt
rm *_snp.vcf.gz
rm *.snp.txt
rm vcf_files_list.txt

# concatenate indel_mutation_table.txt and snp_mutation_table.txt

cat indel_mutation_table.txt snp_mutation_table.txt >> unified_variants.txt
