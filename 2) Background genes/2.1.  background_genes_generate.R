# split into a list
my_genes<-read.delim(file = "my_genes",  header = TRUE, sep=" ")
hg38_genes<-read.delim(file = "hg38_genes",  header = FALSE)

my_genes_list<-split(my_genes,seq(nrow(my_genes)))
names(my_genes_list)<-my_genes$Gene


# false gene generato function
false_gene_generator<-function(genelist,samplesize,genelength,genename){
  # sample the chr and start from hg19 canonical
  # using genelength as seed
  genelength<-as.integer(genelength)
  set.seed(genelength)
  # use sample function, index dataframe and just return chr and start
  false_gene_df<-genelist[sample(nrow(genelist), samplesize), c(1,2)]
  # add the gene length to make the end
  false_gene_df$end<-false_gene_df[,2]+genelength
  # add in gene name, run and unique id
  false_gene_df$genename<-genename
  false_gene_df$runs<-paste("X",c(1:samplesize),sep="")
  false_gene_df$testname<-paste(false_gene_df$runs,false_gene_df$genename, sep="")
  # give appropriate column names
  names(false_gene_df)<-c("Chr","Chrom_start","Chrom_end","genenames","runs","idname")
  
  # return the dataframe
  return(false_gene_df)
}

# loop through
false_gene_list<-lapply(my_genes_list,function(x) false_gene_generator(hg38_genes,1500,x[,2],x[,1]))

# bind up into single data table for intersecting
write.table(do.call(rbind,false_gene_list),file="fake_genes_table", sep="\t", col.names= F, row.names= F, quote= F)
