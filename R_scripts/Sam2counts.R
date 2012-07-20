##R script to obtain counts per transcript when reads have been mapped to cDNAs

##searches working directory for .sam files and summarizes them.

#BEFORE running the lines below change the working directory to 
#the directory with sam files.  If your files end in something
#other than ".sam" change the command below

#get a list of sam files in the working directory
#the "$" denotes the end of line
files <- list.files(pattern="\\.sam$")

#look at files to make sure it is OK
print(files)

#create an empty object to hold our results
results <- NULL

#loop through each file...
for (f in files) {
  print(f) #print the current file
  
  #read the file.  We only care about the third column.
  #also discard the header info (rows starting with "@")
  tmp <- scan(f,what=list(NULL,NULL,""),
            comment.char="@",sep="\t",flush=T)[[3]]
  
  #use table() to count the occurences of each gene.
  #convert to matrix for better formatting later
  tmp.table <- as.data.frame(table(tmp))
  colnames(tmp.table) <- c("gene",f) #get column name specified
  #not needed, in fact a mistake, I think. 
  #tmp.table$gene <- rownames(tmp.table)
  
  #add current results to previous results table, if appropriate
  if (is.null(results)) { #first time through
    results <- as.data.frame(tmp.table) #format
    } else { #not first time through
      results<-merge(results,tmp.table,all=T,
                     by="gene") #combine results
      #rownames(results) <- results$Row.names #reset rownames for next time through
    } #else
  } #for
  rm(list=c("tmp","tmp.table")) #remove objects no longer needed

#summarize mapped and unmapped reads:
print("unmapped")
unmapped <- results[results$gene=="*",-1]
unmapped
results.map <- results[results$gene!="*",]
print("mapped")
mapped <- apply(results.map[-1],2,sum,na.rm=T)
mapped
print("percent mapped")
round(mapped/(mapped+unmapped)*100,1)


write.table(results.map,file="sam2countsResults.tsv",sep="\t",row.names=F)



