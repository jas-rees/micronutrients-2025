#sample branch lengths from Relate
#done so SNPs specified in "calcium_iron_snps"

while read line; do
chr=$(echo $line | awk '{print $2}')
pos=$(echo $line | awk '{print $3}')

while read pop; do 
/home/ssd/jrees/software/relate_v1.1.6_x86_64_static/scripts/SampleBranchLengths/SampleBranchLengths.sh --mu 1.25e-8 -i /home/ssd/jrees/relate/popsize/${pop}/hgdp_poplabels_${chr} -o branch_files/${pop}/hgdp_poplabels_${chr}_${pos} --coal /home/ssd/jrees/relate/popsize/${pop}/hgdp_poplabels.coal --format n --num_samples 200 --first_bp ${pos} --last_bp ${pos} --dist /home/ssd/jrees/relate/popsize/${pop}/hgdp_poplabels_${chr}.dist

done < /home/ssd/jrees/popdata/pops_grouped.txt
done < calcium_iron_snps

#see https://github.com/avaughn271/CLUES2 for full documentation and description of file formats
#generate allele files

#make lists of individuals 
awk '{print $1}' /home/ssd/jrees/relate/prepd_files/chr1_prepd.sample > list_indvs
sed -i '1d' list_indvs
sed -i '1d' list_indvs

while read line; do
echo $line >> list_indvs_both
echo $line >> list_indvs_both
done < list_indvs

#take care below when some group names have "_" - have edited the pops_grouped.txt to avoid issues 
while read line; do
pos=$(echo $line | awk '{print $3}')
chr=$(echo $line | awk '{print $2}')
echo $pos $chr 

while read pop; do
awk -v pop="$pop" '{if ($2 == pop) print $1}' /home/ssd/jrees/relate/hgdp.poplabels > tmp
while read id; do
awk -v id="$id" '{if ($2 == id) print $1}' ${chr}_${pos}_allele_indvs >> ${pop}_${chr}_${pos}_allele_indvs
done < tmp
done < /home/ssd/jrees/popdata/pops_grouped_edit.txt
done < /home/ssd/jrees/clues2_micros/calcium_iron_snps

#check if any alleles are flipped
while read line; do
pos=$(echo $line | awk '{print $3}')
chr=$(echo $line | awk '{print $2}')

status=$(awk -v pos="$pos" '{if ($3 == pos) print $4, $5}' /home/ssd/jrees/relate/prepd_files/${chr}_prepd.haps)
flipped=$(awk -F ';' -v pos="$pos" '{if ($2 == pos) print $8}' /home/ssd/jrees/relate/popsize/Adygei/hgdp_poplabels_${chr}.mut)

echo $gene $chr $pos $status $flipped
echo $gene $chr $pos $status $flipped >> /home/ssd/jrees/clues2_micros/flip_check
done < /home/ssd/jrees/clues2_micros/calcium_iron_snps

#RelatetoCLUES step (see https://github.com/avaughn271/CLUES2)

while read line; do
pos=$(echo $line | awk '{print $3}')
chr=$(echo $line | awk '{print $2}')
while read pop; do
python3 /project/tishkofflab/projects/topmed_alpha/reesjas/software/CLUES2/RelateToCLUES.py --RelateSamples branch_files/${pop}_${chr}_${pos}.newick --DerivedFile allele_files/${pop}_${chr}_${pos}_allele_indvs --out output/${pop}_${chr}_${pos}
done < pops_grouped.txt
done < calcium_iron_snps

#inference step (see https://github.com/avaughn271/CLUES2)
while read line; do
pos=$(echo $line | awk '{print $3}')
chr=$(echo $line | awk '{print $2}')
while read pop; do
for time in 500 1000 1500 2000; do
freq=$(awk -v pos="$pos" -v chr="$chr" -v pop="$pop" '{if ($1 == pop && $2 == chr && $3 == pos) print $4}' population_frequencies)

python3 /project/tishkofflab/projects/topmed_alpha/reesjas/software/CLUES2/inference.py 
--coal /project/tishkofflab/projects/topmed_alpha/reesjas/trpm8_clues/HGDP/data/coal_HGDP/${pop}_hgdp_poplabels.coal
--popFreq ${freq} 
--times /project/tishkofflab/projects/topmed_alpha/reesjas/NEW_clues_data/clues2_micros/output/${pop}_${chr}_${pos}_times.txt 
--out /project/tishkofflab/projects/topmed_alpha/reesjas/NEW_clues_data/clues2_micros/output/${pop}_${chr}_${pos}_${time} 
--tCutoff ${time} --df 450 --noAlleleTraj

done
done < pops_grouped.txt
done < calcium_iron_snps


