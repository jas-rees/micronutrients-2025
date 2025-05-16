#all_micros_summed_pop in format <pop> <sum> where sum is the sum of the pvalues of SNPs with the strongest evidence of positive selection (over all MA-genes)
#below is for FST, analagous for Relate

import pandas as pd
import numpy as np
from scipy.stats import norm

pops=["BantuSouthAfrica_BantuKenya", "Biaka", "Mandenka", "Mbuti", "San", "Bedouin", "Druze", "Mozabite", "Palestinian", "Adygei", "Basque", "BergamoItalian_Tuscan", "French", "Orcadian", "Russian", "Sardinian", "Balochi", "Brahui", "Burusho", "Hazara", "Kalash", "Makrani", "Pathan", "Sindhi", "Uygur", "Dai_Lahu", "Han", "Japanese", "Oroqen_Hezhen_Daur", "Naxi_Yi", "NorthernHan_Tu", "She_Miao_Tujia", "Xibo_Mongolian", "Yakut", "Maya", "Pima", "Surui_Karitiana", "Bougainville", "PapuanHighlands_PapuanSepik"]
micro=["Selenium", "Copper", "Iron", "Magnesium", "Zinc", "Sodium", "Calcium", "Iodine", "Chloride", "Potassium", "Phosphorus", "Manganese", "Molybdenum"]
newfile=[]
for x in pops:
	dist = pd.read_csv("/home/ssd/jrees/significant/fst/sumstat/all/{}_neutral_summed".format(x), header=None)
	mean = np.mean(dist[0])
	std = np.std(dist[0])
	file=pd.read_csv("/home/ssd/jrees/significant/fst/sumstat/all_micros_summed_pop", header=None, sep=" ")
	score = file.loc[file[0]=="{}".format(x)][1].values[0]
	prob = 1- norm(mean, std).cdf(score)
	newfile.append(prob) 
file['Probability'] = newfile
np.savetxt("/home/ssd/jrees/significant/fst/sumstat/all_micros_summed_pop", file, fmt = '%s %f %f',header="Population Sum Probability", comments="")
