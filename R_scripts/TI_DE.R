#analyzing RNAseq for differential expression
#Julin Maloof
#August, 2011
#For Tomato Interns

#the following two lines only need to be run once on your computer
#they install the package that we use for differential expression.
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")


library(edgeR)

#you can see a detailed help file with
vignette("edgeR") #read this later if you want more details

#read in raw count data per gene
counts <- read.delim("sam2countsResults.tsv",row.names=1)

#check the file
head(counts)
summary(counts)
#need to convert NA to 0 counts
counts[is.na(counts)] <- 0

#convert data to a form that edgeR wants
data <- DGEList(counts=counts,group=rep(c("MT","WT"),times=c(3,4)))
  
data$samples
                
#normalize library
data <- calcNormFactors(data)

data$samples
                
#estimate overdispersion
#important so that the correct model is fit
data <- estimateCommonDisp(data)

#calculate DE genes
DEtest <- exactTest(data,pair=c("WT","MT"))
head(DEtest)

#create a table of the results, with multiple testing correction
results <- topTags(DEtest,n=Inf)
                
head(results$table,n=30) #top 30 genes

#how many genes are DE?

sum(results$table$FDR<.01) 

#how many genes in each direction?
summary(decideTestsDGE(DEtest,p.value=.01))


#plot the results
#first create a table of DE  to highlight those with p < 0.01
sig.genes <- rownames(results$table[results$table$FDR<0.01,])

head(sig.genes)

plotSmear(data,de.tags=sig.genes)

#what are the genes that are misexpressed?
#for this we need to add some annotation

#need to get rid of 2nd part of "GSVIVT01017617001|pacid:17827756"
row.names(results$table) <- sub("\\|.*","",row.names(results$table))

annotation <- read.delim("Vvinifera_145_annotation_info.txt", head=F)
head(annotation)
head(results$table)

results.annotated <- merge(results$table,annotation,by.x="row.names",by.y="V1",all.x=T,sort=F)

head(results.annotated,n=30)

#perhaps easier to view in excel
write.table(results.annotated,"WT_MT_DE.tsv",sep="\t",row.names=F)




