#filter vcfs for SHAPEIT2
#SLURM_ARRAY_TASK_ID is what is used in this script to indicate chromosome number

/home/ssd/jrees/software/vcftools/bin/vcftools 
--gzvcf /home/ssd/jrees/rawdata/vcfs/hgdp_wgs.20190516.full.chr${SLURM_ARRAY_TASK_ID}.vcf.gz 
--bed /home/ssd/jrees/rawdata/masks/mask_chr${SLURM_ARRAY_TASK_ID}.bed 
--min-alleles 2 
--max-alleles 2 
--max-missing-count 185 
--recode 
--out chr${SLURM_ARRAY_TASK_ID}_preshapeit

gzip chr${SLURM_ARRAY_TASK_ID}_preshapeit.recode.vcf

#run SHAPEIT2
/home/ssd/jrees/software/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit 
--input-vcf /home/ssd/jrees/shapeit/preshapeit/chr${SLURM_ARRAY_TASK_ID}_preshapeit.recode.vcf.gz \
-M /home/ssd/jrees/rawdata/maps/map_chr${SLURM_ARRAY_TASK_ID}.txt \
-O chr${SLURM_ARRAY_TASK_ID} \
-T 88 \
--window 0.3 \
--states 200 \
--output-log chr${SLURM_ARRAY_TASK_ID}.log


# generate hg38 mask in fa format
#SLURM_ARRAY_TASK_ID is what is used in this script to indicate chromosome number

tar -zxvf /home/ssd/jrees/rawdata/ancestrals/homo_sapiens_ancestor_GRCh38.tar.gz homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${SLURM_ARRAY_TASK_ID}.fa 
cd homo_sapiens_ancestor_GRCh38

##Turn every character in ancestral that isn't 'X' into a 'N' for each chromosome
##gives a single line file of Ns

sed -i '1d' homo_sapiens_ancestor_${SLURM_ARRAY_TASK_ID}.fa 
tr -d '\n' < homo_sapiens_ancestor_${SLURM_ARRAY_TASK_ID}.fa > tmp${SLURM_ARRAY_TASK_ID} && rm homo_sapiens_ancestor_${SLURM_ARRAY_TASK_ID}.fa
tr -c 'X' 'N' < tmp${SLURM_ARRAY_TASK_ID} > chr${SLURM_ARRAY_TASK_ID}.fa && rm tmp${SLURM_ARRAY_TASK_ID} 
sed -i -e '$a\' chr${SLURM_ARRAY_TASK_ID}.fa

##replaces regions with Ps as specified in the bed mask file

k=$(cat /home/ssd/jrees/rawdata/masks/mask_chr${SLURM_ARRAY_TASK_ID}.bed | wc -l)

for x in `eval echo {2..$k}`; do

a=$(awk -v x="$x" 'NR==x{print $2}' /home/ssd/jrees/rawdata/masks/mask_chr${SLURM_ARRAY_TASK_ID}.bed);
b=$(awk -v x="$x" 'NR==x{print $3}' /home/ssd/jrees/rawdata/masks/mask_chr${SLURM_ARRAY_TASK_ID}.bed);

awk -v a="$a" -v b="$b" '{t=substr($0,a,b-a+1); gsub(/./,"P",t); print substr($0,1,a-1) t substr($0,b+1)}' chr${SLURM_ARRAY_TASK_ID}.fa > tmp${SLURM_ARRAY_TASK_ID} && mv tmp${SLURM_ARRAY_TASK_ID} chr${SLURM_ARRAY_TASK_ID}.fa; 
done

##cuts file at 100 each line

fold -w100 chr${SLURM_ARRAY_TASK_ID}.fa > mask_chr${SLURM_ARRAY_TASK_ID}.txt
rm chr${SLURM_ARRAY_TASK_ID}.fa

##gives header to fasta file 

n=$(grep -o "N" mask_chr${SLURM_ARRAY_TASK_ID}.txt | wc -l)
p=$(grep -o "P" mask_chr${SLURM_ARRAY_TASK_ID}.txt | wc -l)

sed  -i "1i >21 dna:chromosome chromosome: GRCh38:21 [$(expr $p + $n) bases - $n are N, $p are P]" mask_chr${SLURM_ARRAY_TASK_ID}.txt
mv mask_chr${SLURM_ARRAY_TASK_ID}.txt /home/ssd/jrees/relate/dist_prep/mask_chr${SLURM_ARRAY_TASK_ID}.fa

# prepare files
#SLURM_ARRAY_TASK_ID is what is used in this script to indicate chromosome number
#see https://myersgroup.github.io/relate/ for full documentation

/home/ssd/jrees/software/relate_v1.1.2_x86_64_static/scripts/PrepareInputFiles/PrepareInputFiles.sh 
--haps /home/ssd/jrees/shapeit/chr${SLURM_ARRAY_TASK_ID}.haps \
--sample /home/ssd/jrees/shapeit/chr${SLURM_ARRAY_TASK_ID}.sample \
--ancestor /home/ssd/jrees/rawdata/ancestrals/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${SLURM_ARRAY_TASK_ID}.fa \
--mask /home/ssd/jrees/relate/dist_prep/mask_chr${SLURM_ARRAY_TASK_ID}.fa \
--poplabels /home/ssd/jrees/rawdata/hgdp.poplabels \
-o chr${SLURM_ARRAY_TASK_ID}_prepd

#replace paths to files appropriately
#multithreaded

/home/ssd/jrees/software/relate_v1.1.6_x86_64_static/scripts/RelateParallel/RelateParallel.sh \
--mode All \
-m 1.25e-8 \
-N 30000 \
--haps /home/ssd/jrees/relate/prepd_files/chr${SLURM_ARRAY_TASK_ID}_prepd.haps.gz \
--sample /home/ssd/jrees/relate/prepd_files/chr${SLURM_ARRAY_TASK_ID}_prepd.sample.gz \  
--map /home/ssd/jrees/rawdata/maps/map_chr${SLURM_ARRAY_TASK_ID}.txt \
--annot /home/ssd/jrees/relate/prepd_files/chr${SLURM_ARRAY_TASK_ID}_prepd.annot \
--dist /home/ssd/jrees/relate/prepd_files/chr${SLURM_ARRAY_TASK_ID}_prepd.dist.gz \
--seed 1 \
-o chr${SLURM_ARRAY_TASK_ID} \
--threads 18 \ 

#replace pops_grouped.txt file with appropriate path to list of populations
#replace paths to input (haps, sample) files appropriately
#multithreaded

while read pop; do 
for chr in {1..22}; do
/home/ssd/jrees/software/relate_v1.1.6_x86_64_static/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \ 
-i /home/ssd/jrees/relate/mut_anc/chr${SLURM_ARRAY_TASK_ID} \
-m 1.25e-8 \
--poplabels hgdp.poplabels
--seed 1 \
--threshold 0 \
--pop_of_interest ${pop} \ 
--threads 18 \
-o hgdp_poplabels \
done  
done < pops_grouped.txt

#replace path to list of populations
#replace path to input files (mut, anc)

while read pop; do
for chr in {1..22}; do

mkdir ${pop}

/home/ssd/jrees/software/relate_v1.1.6_x86_64_static/scripts/DetectSelection/DetectSelection.sh \
-i /home/ssd/jrees/relate/popsize/${pop}/hgdp_poplabels_chr${chr} \ 
-o ${pop}/${pop}_chr${chr}_selection \
-m 1.25e-08 \
-years_per_gen 28 \ 

done 
done < /home/ssd/jrees/relate/pops
