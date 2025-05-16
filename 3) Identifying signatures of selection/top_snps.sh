#SNPs with the strongest evidence of selection in each MA-gene
# below for Relate
# relate/micro_per_pop/${pop}/chr${chrom}_micronutrients in format <micronutrient> <gene> <position> <rsID> <relate_score>

while read micro; do
while read line; do
while read pop; do

gene=$(echo ${line} | awk '{print $2}')
chrom=$(echo ${line} | awk '{print $3}')

top_sig=$(grep ${gene} /home/ssd/jrees/significant/relate/micro_per_pop/${pop}/chr${chrom}_micronutrients | sort -k 5 -n | head -n 1 | awk '{print $5}')

echo ${pop} ${micro} ${gene} ${top_sig} >> /home/ssd/jrees/significant/relate/micro_per_pop/top_snps_per_pop

done < /home/ssd/jrees/popdata/pops_grouped.txt
done < /home/ssd/jrees/micronutrients/cut_by_micro/${micro}_all
done < /home/ssd/jrees/micronutrients/micronutrients_list

#below for Fst
# fst/micro_per_pop/${pop}/chr${chrom}_micronutrients in format <micronutrient> <gene> <position> <fst_score>

while read micro; do
while read line; do
while read pop; do

gene=$(echo ${line} | awk '{print $2}')
chrom=$(echo ${line} | awk '{print $3}')

top_sig=$(grep ${gene} /home/ssd/jrees/significant/fst/micro_per_pop/${pop}/chr${chrom}_micronutrients | sort -k 4 -r | head -n 1 | awk '{print $4}')

echo ${pop} ${micro} ${gene} ${top_sig} >> /home/ssd/jrees/significant/fst/micro_per_pop/top_snps_per_pop

done < /home/ssd/jrees/popdata/pops_grouped.txt
done < /home/ssd/jrees/micronutrients/cut_by_micro/${micro}_all
done < /home/ssd/jrees/micronutrients/micronutrients_list
