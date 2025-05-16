#selecting top snps (example given for relate)
#analagous done for neutral snps

for pop in africa europe asia america; do
for year in 1 5 10 40; do
while read line; do
seed=$(echo ${line} | awk '{print $1}')
awk -v seed="$seed" '{if ($1 == seed) print $0}' ${pop}_${year}_all_sites_relate_prob_ecdf | sort -k 8 -g | head -n 1 | awk '{print $0}' >> top_snps/${pop}_${year}_relate_top_snps
done < ../${year}kya/${pop}/all_success.txt
done
done

#and summing 
#analagous done for neutral snps

for pop in africa europe asia america; do
for year in 1 5 10 40; do
for n in 10 20 40 60; do 
for i in {1..1000}; do 
awk '{print $8}' top_snps/${pop}_${year}_relate_top_snps | sed '1d' | shuf -n ${n} | awk '{sum += $1; } END {print sum;}' >> sums/${pop}_${year}_relate_${n}_top_snps
done
done
done
done

#generating gene sets of differing proportion of SNPs under selection
for year in 1 5 10 40; do
for pop in africa europe asia america; do
for n in 10 20 40 60; do
for perc in 0.2 0.4 0.6 0.8; do
x=$(expr 1-$perc | bc)
s=$(expr $n*$perc | bc)
ns=$(expr $n*$x | bc)
for i in {1..1000}; do
selected=$(awk '{print $8}' top_snps/${pop}_${year}_relate_top_snps | shuf -n ${s} | awk '{ sum+= $1; } END {print sum;}')
neutral=$(awk '{print $8}' neutral/${pop}_${year}_relate_top_snps | sed '1d' | shuf -n ${ns} | awk '{ sum+= $1; } END {print sum;}')
sum=$(expr $selected+$neutral | bc)
echo ${sum} >> percentage/sum_${n}_${perc}_${year}kya_${pop}_relate
done
done
done
done
done

