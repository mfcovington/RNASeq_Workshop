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
DEtest <- exactTest(data,pair=c("Heinz","M82"))
head(DEtest)

#create a table of the results, with multiple testing correction
results <- topTags(DEtest,n=Inf)
                
head(results,n=30) #top 30 genes

#how many genes are DE?

sum(results$table$adj.P.Val<.01) 

#how many genes in each direction?
summary(decideTestsDGE(DEtest,p.value=.01))


#plot the results
#first create a table of DE  to highlight those with p < 0.01
sig.genes <- rownames(results$table[results$table$adj.P.Val<0.01,])

head(sig.genes)

plotSmear(data,de.tags=sig.genes)

#what are the genes that are misexpressed?
#for this we need to add some annotation
annotation <- read.delim("ITAG2.3_all_Arabidopsis_annotated.tsv")
head(annotation)
head(results)

results.annotated <- merge(results$table,annotation,by.x="row.names",by.y="ITAG",all.x=T,sort=F)

head(results.annotated,n=30)

#perhaps easier to view in excel
write.table(results.annotated,"M82_HZ_DE.tsv",sep="\t",row.names=F)




