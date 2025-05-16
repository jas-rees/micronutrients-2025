#!/bin/bash
# background gene tables for each micronutrients
# micronutrients_list is format <micronutrient>
# ${micro}_all is format <gene>

while read micro; do
while read line; do
gene=$(echo ${line} | awk '{print $2}')
grep -w ${gene} fake_genes_table > ${micro}_fake_genes_table
done < ${micro}_all
done < micronutrients_list

# extract values of Relate and Fst for all background genes
# pops_grouped.txt is in format <population>
# below is for Relate (uses .sele output from Relate for each population)

while read micro; do
while read pop; do
while read line; do
gene=$(echo ${line} | awk '{print $4}')
num=$(echo ${line} | awk '{print $5}')
chr=$(echo ${line} | awk '{print $1}')

lower=$(echo ${line} | awk '{print $2}')
upper=$(echo ${line} | awk '{print $3}')
downstream=$((${lower} - 10000))
upstream=$((${upper} + 10000))

awk -v gene="$gene" -v num="$num" -v downstream="$downstream" -v upstream="$upstream" '{if ($1 > downstream && $1 < upstream) print gene,num,$1,$2,$35}' relate/selection/${pop}/${pop}_${chr}_selection.sele >> ${pop}/${micro}_fake_genes_relate

done < ${micro}_fake_genes_table
done < pops_grouped.txt
done < micronutrients_list

# below is for Fst

while read micro; do
while read pop; do
while read line; do
gene=$(echo ${line} | awk '{print $4}')
num=$(echo ${line} | awk '{print $5}')
chr=$(echo ${line} | awk '{print $1}' | cut -c 4,5)

lower=$(echo ${line} | awk '{print $2}')
upper=$(echo ${line} | awk '{print $3}')
downstream=$((${lower} - 10000))
upstream=$((${upper} + 10000))

awk -v gene="$gene" -v num="$num" -v downstream="$downstream" -v upstream="$upstream" '{if ($2 > downstream && $2 < upstream) print gene,num,$2,$3}' fst/${pop}*Yoruba/${pop}*Yoruba_${chr}_masked >> ${micro}_fake_genes_fst

done < ${micro}_fake_genes_table
done < pops_grouped.txt
done < micronutrients_list
