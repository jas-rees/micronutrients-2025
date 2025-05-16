# for all micronutrient genes ("all_genes"), takes the strongest evidence of positive selection in a random background-gene (corresponding to the MA-gene) and submits to tmp
# sums total and submits to all/${pop}_neutral_summed
# repeats 1000 times
# for loop for each population
# below for fst, analagous for Relate

while read pop; do 
for i in {1..1000}; do
while read gene; do
grep -w ${gene} /home/ssd/jrees/fake_genes/populations/${pop}/*_fake_genes_fst_top_snps | shuf -n 1 | awk '{print $4}' >> tmp
done < all_genes
awk '{s+=$1} END {print s}' tmp >> /home/ssd/jrees/significant/fst/sumstat/all/${pop}_neutral_summed
rm tmp
done
done < /home/ssd/jrees/popdata/pops_grouped.txt


# for each MA-gene in a micronutrient category (${micro}_all), takes the strongest evidence of positive selection in a random background-gene (corresponding to the MA-gene) and submits to$
# sums total and submits to ${micro}/${pop}_neutral_summed
# repeats 1000 times                                                                 
# for loop for each population   
# for loop for each micronutrient
# below for fst, analagous for Relate

while read micro; do
while read pop; do 
for i in {1..1000}; do
while read line; do
gene=$(echo $line | awk '{print $2}')
grep -w ${gene} /home/ssd/jrees/fake_genes/populations/${pop}/*_fake_genes_fst_top_snps | shuf -n 1 | awk '{print $4}' >> tmp
done < ${micro}_all
awk '{s+=$1} END {print s}' tmp >> /home/ssd/jrees/significant/fst/sumstat/${micro}/${pop}_neutral_summed
rm tmp
done
done < /home/ssd/jrees/popdata/pops_grouped.txt
done < /home/ssd/jrees/micronutrients/micronutrients_list

