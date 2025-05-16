# used to evaluate where MA-genes fall in the distribution of their corresponding background genes
# snp_counts_micronutrient in the format <Micronutrient> <Gene> <SNP count> (SNP count inferred from Yoruba)
# ${genes} files list SNP counts of all background-genes for the corresponding MA-gene

import pandas as pd
import numpy as np
from scipy.stats import norm

newfile=[]
file = pd.read_csv("np_counts_micronutrient", sep=' ')
rows=len(file)

for i in range(rows):
	gene=file["Gene"].values[i]
	score=file["SNP_count"].values[i]
	dist = pd.read_csv("genes/{}".format(gene), header=None)
	mean = np.mean(dist[0])
	std = np.std(dist[0])
	prob = norm(mean, std).cdf(score)
	newfile.append(prob)

file['CDF']= newfile
np.savetxt("/home/ssd/jrees/fake_genes/CDF", file, fmt= '%d %d %d %f', header="Micronutrient Gene SNP CDF", comments="")
