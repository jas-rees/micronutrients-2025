#filter vcfs for biallelic sites, removing indels
for chr in {1..22}; do
/home/ssd/jrees/software/vcftools/bin/vcftools --gzvcf /home/ssd/jrees/raw_vcfs/hgdp_wgs.20190516.full.chr${chr}.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2  --recode --recode-INFO-all --stdout | gzip -c >  biallelic_chr${chr}.vcf.gz 
done

#run over population combinations (given by combos.txt)
while read line; do
pop1=$(echo $line | cut -f 1 -d "*")
pop2=$(echo $line | cut -f 2 -d "*")
for chr in {1..22}; do
/home/ssd/jrees/software/vcftools/bin/vcftools --gzvcf biallelic_chr${chr}.vcf.gz --weir-fst-pop /home/ssd/jrees/popdata/${pop1}_list.txt --weir-fst-pop /home/ssd/jrees/popdata/${pop2}_list.txt --out ${line}/${line}_${chr}
done
done < combos.txt

# To remove -nan values and turn negative sites to 0 and remove e notation
while read pop; do
sed '/-nan/d'  /home/ssd/jrees/fst/${pop}*Yoruba/${pop}*Yoruba_2.weir.fst | awk '{$3=($3<0)?0:$3}1' | awk '{$3=sprintf("%.24f",$3)}7' > /home/ssd/jrees/fst/${pop}*Yoruba/${pop}*Yoruba_2_edited.weir.fst
# remove header
sed -i -e "1d" /home/ssd/jrees/fst/${pop}*Yoruba/${pop}*Yoruba_2_edited.weir.fst
done < /home/ssd/jrees/popdata/pops_grouped.txt

#Apply masking
for chr in {1..22}; do
while read pop; do
while read line; do
start=$(echo $line | awk '{print $2}')
end=$(echo $line | awk '{print $3}')
awk -v start=$start -v end=$end '{if ($2 >= start && $2 <= end) print $0}' /home/ssd/jrees/fst/${pop}*Yoruba/${pop}*Yoruba_${chr}_edited.weir.fst  >> /home/ssd/jrees/fst/${pop}*Yoruba/${pop}*Yoruba_${chr}_masked.weir.fst
done < /home/ssd/jrees/mask/chr${chr}_mask
done < /home/ssd/jrees/popdata/pops_grouped.txt
done
