# calculate empirical p-values for the highest ranking SNP in each MA-gene
# dist specifies the empirical distribution of selection values 
# below for Fst


emp_prob<-function(x,num){
+ sum(num<=x)/length(x)}

top<-read.csv("/home/ssd/jrees/significant/fst/micronutrients/top_snps_each_pop", header=TRUE, sep=" ")

for (p in c("San", "BantuSouthAfrica_BantuKenya", "Mbuti", "Biaka", "Mandenka", "Mozabite", "Palestinian", "Druze", "Bedouin", "Adygei", "BergamoItalian_Tuscan", "Sardinian", "Basque", "French", "Orcadian", "Russian", "Makrani", "Sindhi", "Balochi", "Brahui", "Hazara", "Pathan", "Burusho","Kalash", "Uygur", "Xibo_Mongolian", "Oroqen_Hezhen_Daur", "Yakut", "Japanese", "Han", "NorthernHan_Tu", "She_Miao_Tujia", "Naxi_Yi", "Dai_Lahu", "Pima", "Maya", "Surui_Karitiana", "PapuanHighlands_PapuanSepik", "Bougainville")){
dist<-readLines(paste0(p,"*Yoruba/all"))
for (i in 1:277) {
top[i,p] <- emp_prob(dist,top[i,p])
}
}
write.table(top, file = "top_snps_p_values_fst", sep= "\t", col.names = FALSE)

#below fo Relate
emp_prob<-function(x,num){
+ sum(num<=x)/length(x)}

top<-read.csv("/home/ssd/jrees/significant/relate/micronutrients/top_snps_each_pop", header=TRUE, sep=" ")

for (p in c("San", "BantuSouthAfrica_BantuKenya", "Mbuti", "Biaka", "Yoruba", "Mandenka", "Mozabite", "Palestinian", "Druze", "Bedouin", "Adygei", "BergamoItalian_Tuscan", "Sardinian", "Basque", "French", "Orcadian", "Russian", "Makrani", "Sindhi", "Balochi", "Brahui", "Hazara", "Pathan", "Burusho","Kalash", "Uygur", "Xibo_Mongolian", "Oroqen_Hezhen_Daur", "Yakut", "Japanese", "Han", "NorthernHan_Tu", "She_Miao_Tujia", "Naxi_Yi", "Dai_Lahu", "Pima", "Maya", "Surui_Karitiana", "PapuanHighlands_PapuanSepik", "Bougainville")){
dist<-readLines(paste0(p,"_pvalues"))
for (i in 1:277) {
top[i,p] <- emp_prob(dist,top[i,p])
}
}
write.table(top, file = "top_snps_p_values_relate", sep= "\t", col.names = FALSE)
