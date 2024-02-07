#### Convert Qiime/DADA2 output into tables to be used in R ####

#### STEP 1: Load Libraries ####
library(plyr) # Load plyr library into R environment

#### STEP 2: Import Data ####

count<- read.delim("/Users/briannepalmer/Downloads/featuretables/feature-table_dada2.tsv", header=FALSE, row.names=NULL, stringsAsFactors=FALSE, na.strings = "n/a")

#### STEP 3: Reclassify and Reformat ####
colnames(count)<-as.character(unlist(count[2,])) # makes column names the values from the second row 
count1<-count[-c(1,2),] # removes the first 2 rows 
colnames(count1)[1]<-"Feature.ID" # change the name of column 1
head(count1[1:2,]) # check dataframe structure, check that the row and column names are correct 
x<-dim(count1)[2];x # number of columns (should be total sample number + 1)

# Convert the dataframe so the data is "numeric" instead of a "character"
count2<-count1
x<-dim(count2)[2];x # number of columns
count2[2:x] <- lapply(count2[2:x], function(x) as.numeric(as.character(x)))

##### STEP 4: Get base stats and generate a table called "dataset_info" ####

seq_total<-apply(count2[2:x],2,sum) #number of sequences per sample
OTU_count<-colSums(count2[2:x]>0) # ASVs per sample
OTU_single<-colSums(count2[2:x]==1) # ASVs with only 1 seq -- these are usually not reliable 
OTU_double<-colSums(count2[2:x]==2) # ASVs that have only 2 seqs
OTU_true<-colSums(count2[2:x]>2) # Number of ASVs with >2 seqs

# Option to remove global singletons (this is good to do!)
rowsum<-apply(count2[2:x],1,sum) # Number of sequences per ASV (for whole dataset)
count.no1 = count2[ rowsum>1, ] # Remove rowsums =< 1 (only one sequence for 1 ASV)
dim(count2)[1] - dim(count.no1)[1] # Total number of global singletons removed from data

dataset_info1 <- data.frame(seq_total,OTU_count,OTU_single,OTU_double,OTU_true)

write.table(dataset_info1, file= "path/to/where/you/keep/output", quote=FALSE, sep="\t", row.names=FALSE)

# Totals for whole dataset:
sum(seq_total) # total number of sequences over whole dataset
dim(count2)[1] # total unique OTUs over whole dataset

##### STEP 5: Import and Merge Taxonomy 
# Import taxonomy file and join with count data:
taxname<-read.delim("path/to/taxonomy/file", header=TRUE, row.names=NULL) # taxonomy file 
counts_wtax<-join(count.no1, taxname, by="Feature.ID", type="left", match="first") # merge with ASV file created above 

# Write output table with ASV cluster information and taxonomy IDs:
write.table(counts_wtax, file= "path/to/table/output", quote=FALSE, sep="\t", row.names=FALSE)


