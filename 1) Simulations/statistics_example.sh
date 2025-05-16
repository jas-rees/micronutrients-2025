#iHS and nSL filtering example

while read line; do
/home/ucsaree/software/vcftools/bin/vcftools --gzvcf ../1kya/africa/${line}.vcf.gz --min-alleles 2 --max-alleles 2 --maf 0.05 --max-maf 0.95 --recode --out ../1kya/africa/${line}_preihs
done < ../1kya/africa/seed_africa.txt

#iHS and nSL calculation example
while read line; do
/home/ucsaree/software/selscan-linux-1.3.0/selscan --ihs --ihs-detail --vcf ../1kya/africa/${line}_preihs.recode.vcf.gz --map ../1kya/africa/maps/ihs_map_africa_${line}.map --trunc-ok --out ../1kya/africa/${line}_iHS
done < ../1kya/africa/seed_africa.txt

while read line; do
/home/ucsaree/software/selscan-linux-1.3.0/selscan --nsl --vcf ../1kya/africa/${line}_preihs.recode.vcf.gz --trunc-ok --out ../1kya/africa/${line}_nSL
done < ../1kya/africa/seed_africa.txt

#iHS and nSL normalisation example
/home/ucsaree/software/selscan-linux-1.3.0/norm --ihs --files ../1kya/africa/${line}_iHS*
/home/ucsaree/software/selscan-linux-1.3.0/norm --nsl --files ../1kya/africa/${line}_nSL*

#XP-EHH and XP-nSL pre-processing examples
#merging population vcfs so the same sites are present in the vcfs when we extract

/home/ucsaree/software/htslib-1.11/bgzip ../1kya/africa/${line}.vcf.gz
/home/ucsaree/software/bcftools-1.11/bcftools index ../1kya/africa/${line}.vcf.gz

while read line; do
/home/ucsaree/software/bcftools-1.11/bcftools merge ../1kya/africa/${line}.vcf.gz ../neutral/europe/${line}.vcf.gz ../neutral/asia/${line}.vcf.gz ../neutral/america/${line}.vcf.gz -O z -0 --force-samples -o merged_${line}.vcf.gz
done < ../1kya/africa/seed_africa.txt

for x in africa eastasia europe america; do
while read seed; do
/home/ucsaree/software/vcftools/bin/vcftools --gzvcf merged_${line}.vcf.gz /home/ucsaree/Scratch/slim/${x}_indvs.txt --min-alleles 2 --max-alleles 2 --recode --out africa_${seed}_prexpehh
done < ../1kya/africa/seed_africa.txt
done

#XP-EHH and XP-nSL calculation example
while read seed; do
/home/ucsaree/software/selscan-linux-1.3.0/selscan --xpehh --vcf africa_${seed}_prexpehh.recode.vcf.gz --vcf-ref europe_${seed}_prexpehh.recode.vcf.gz --map ../1kya/africa/maps/xpehh_map_africa_${line}.map --trunc-ok --out AFR_EUR_xpehh
done < ../1kya/africa/seed_africa.txt

while read seed; do
/home/ucsaree/software/selscan-linux-1.3.0/selscan --xpnsl --vcf africa_${seed}_prexpehh.recode.vcf.gz --vcf-ref europe_${seed}_prexpehh.recode.vcf.gz --trunc-ok --out AFR_EUR_xpnsl
done < ../1kya/africa/seed_africa.txt

#XP-EHH and XP-nSL normalisation example
/home/ucsaree/software/selscan-linux-1.3.0/norm --xpehh --files AFR_EUR_xpehh*
/home/ucsaree/software/selscan-linux-1.3.0/norm --xpnsl --files AFR_EUR_xpnsl

#Fst example
while read seed; do
/home/ucsaree/software/vcftools/bin/vcftools --gzvcf merged_${seed}.vcf.gz --weir-fst-pop africa_indvs.txt --weir-fst-pop europe_indvs.txt --out AFR_EUR_${seed}
done < ../1kya/africa/seed_africa.txt

#Relate example
while read seed; do
/home/ucsaree/software/shapeit/bin/shapeit --input-vcf ../1kya/africa/africa_${seed}_prexpehh -M ../recombination/maps/map1_${seed}.txt -O haps_sample/${seed} -T 18 --window 0.3
/home/ucsaree/software/relate_v1.0.17_x86_64_static/scripts/RelateParallel/RelateParallel.sh --mode All -m 1.25e-8 -N 30000 --haps haps_sample/${seed}.haps --sample haps_sample/${seed}.sample --map ../recombination/maps/map1_${line}.txt -o ${line} --threads 2
/home/ucsaree/software/relate_v1.0.17_x86_64_static/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i ${seed} --poplabels poplabels.txt --threshold 0 --threads 2 --pop_of_interest AFR -o popsize/${seed}_popsize
/home/ucsaree/software/relate_v1.0.17_x86_64_static/scripts/DetectSelection/DetectSelection.sh -i popsize/${seed} -o selection/AFR_${seed} -m 1.25e-8 --years_per gen 28
done < ../1kya/africa/seed_africa.txt


