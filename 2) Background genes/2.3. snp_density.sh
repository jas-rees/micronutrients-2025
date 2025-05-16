#!/bin/bash

# extract and filter Yoruba vcfs
for chr in {1..22}; do
/home/ssd/jrees/software/vcftools/bin/vcftools --gzvcf hgdp_wgs.20190516.full.chr${chr}.vcf.gz --bed chr${chr}_mask --min-alleles 2 --max-alleles 2 --keep popdata/Yoruba_list.txt --recode --out Yoruba_vcfs/chr${chr}_filtered_Yoruba
done

#SNP count for each MA-gene
while read line; do
micro=$(echo ${line} | awk '{print $1}')
gene=$(echo ${line} | awk '{print $2}')
chr=$(echo ${line} | awk '{print $3}')

lower=$(echo ${line} | awk '{print $4}')
upper=$(echo ${line} | awk '{print $5}')
downstream=$((${lower} - 10000))
upstream=$((${upper} + 10000))

num=$(awk -v downstream="$downstream" -v upstream="$upstream" '{if ($2 > downstream && $2 < upstream) print $0}' Yoruba_vcfs/chr${chr}_filtered_Yoruba.recode.vcf | wc -l)
echo $micro $gene $num >> snp_counts_micronutrient
done < all_micronutrients

#snp count for all background genes generated for each MA-gene
while read micro; do
while read line; do
gene=$(echo ${line} | awk '{print $4}')
gene_num=$(echo ${line} | awk '{print $5}')
chr=$(echo ${line} | awk '{print $1}')

lower=$(echo ${line} | awk '{print $2}')
upper=$(echo ${line} | awk '{print $3}')
downstream=$((${lower} - 10000))
upstream=$((${upper} + 10000))

num=$(awk -v downstream="$downstream" -v upstream="$upstream" '{if ($2 > downstream && $2 < upstream) print $0}' Yoruba_vcfs/${chr}_filtered_Yoruba.recode.vcf | wc -l)
echo $gene $gene_num $num >> snp_counts/snp_counts_fake_genes_${micro}
done < ${micro}_fake_genes_table
done < micronutrient_list

#remove any background genes with 0 SNPS
while read micro; do 
awk '{if ($3 != "0") print $0}' snp_counts_fake_genes_${micro} > snp_counts_fake_genes_${micro}_filtered
done < micronutrients_list

#take the background genes with the closest matching SNP density to their respective background genes
while read line; do
micro=$(echo ${line} | awk '{print $1}')
gene=$(echo ${line} | awk '{print $2}')

length=$(grep -w ${gene} snp_counts_fake_genes_${micro}_filtered | wc -l)

grep -w ${gene} snp_counts_micronutrient >> snp_counts_fake_genes_${micro}_filtered 
true_pos=$(grep -w ${gene} snp_counts_fake_genes_${micro}_filtered | sort -k 3 -n | grep -n ${micro} | cut -f 1 -d ":")

echo $true_pos $length > tmp 

do=$(awk  '{if ($1 < $2/2) print "head -n 1000"; else if ($1 > $2/2) print "tail -n 1000" }' tmp)

sed -i "/${micro}/d" snp_counts_fake_genes_${micro}_filtered

echo "grep -w ${gene} snp_counts_fake_genes_${micro}_filtered | sort -k 3 -n | ${do} > snp_counts_fake_genes_${micro}_keep" >> instructions.sh

echo ${gene} ${micro} ${true_pos} ${length} >> sanity_check

done < all_micronutrients_genes 

#submit instructions file
bsub < instructions.sh

#filter to only include background genes as given in snp_counts_fake_genes_${micro}_keep
while read micro; do
while read line; do
gene=$(echo ${line} | awk '{print $1}')
num=$(echo ${line} | awk '{print $2}')
grep -w "$gene" ${micro}_fake_genes_table | grep -w "$num" > ${micro}_fake_genes_table_keep
done < snp_counts/snp_counts_fake_genes_${micro}_keep
done < micronutrients_list

#finally, get the distribution of SNP counts of all retained background genes
while read micro; do
while read line; do

gene=$(echo ${line} | awk '{print $2}')

grep -w ${gene} snp_counts/snp_counts_fake_genes_${micro}_keep | awk '{print $3}' > genes/${gene}
done < micro_all
done < micronutrients_list
