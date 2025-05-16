#given a file like "all_sites_relate" which gives the calculated stats value (here, it is the Relate statistic) for all simulated SNPs and the neutral distribution of the same statistic in the same population (as calculated from the neutral simulations), we can calculate the position of each SNP in the tail
#then can identify which of the truly selected SNPs fall within e.g., 5% tail 

import pandas as pd
import numpy as np
from scipy.stats import norm
from statsmodels.distributions.empirical_distribution import ECDF

pops=["africa", "europe", "asia", "america"]
for x in pops:
	newfile=[]
	dist = pd.read_csv("neutral/distributions/{}_40_relate".format(x), header=None)
	ecdf = ECDF(dist[0])
	file = pd.read_csv("../monogenic_power/1kya/{}/all_sites_relate".format(x), sep=' ')
	rows = len(file)
	for i in range(rows):
		score = file["log10Pvalue"].values[i]
		prob = ecdf(score)
		newfile.append(prob) 
	file['Probability'] = newfile
	np.savetxt("{}_40_all_sites_relate_prob_ecdf".format(x), file, fmt = '%d %d %f %f %d %d %f %f',header="Seed SelectedSite Freq Selection Position Label log10PValue Probability", comments="")
