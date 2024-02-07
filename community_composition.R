# microbial community composition 

# this starts the same say as alpha and beta diversity, but here we will create plots to look at the overall community composition 

# load libraries 
library(tidyverse) # for general coding 
library(phyloseq) # using phyloseqs objects is easy to manipulate large datasets
library(microeco) # a package full of microbial ecology tools 
library(file2meco) # a package to change phyloseq to microeco objects 
library(ampvis2) # used to make heatmaps and other general calulations, I like how easy it is to calcualte alpha diversity 
library(ggpubr) # needed for graph manipulation 
library(rstatix) # needed for bar plot "add significance"


# upload data 

data <- read.csv("path/to/data")

otu.table <- data[,c(2:55)] # change this to only include the ASV (otu) data
otu.table <- as.matrix(otu.table) # make matrix to change rownames 
rownames(otu.table) <- data$Feature.ID # make FeatureID the rownames 
otu.table <- data.frame(otu.table) # make dataframe again
otu.table <- otu_table(otu.table, taxa_are_rows = TRUE) # converts into phyloseq OTU table 

metadata <- read.csv("path/to/metadata", sep = ",") # metadata file 
rownames(metadata) <- metadata$SampleID # make SampleID rownames 
sam.data <- as.data.frame(metadata) # make dataframe
sam.data.ps <- sample_data(sam.data) # phyloseq object 

tax.table <- data[,c(56:63)] # change this to only include the taxonomy 
tax.table <- as.matrix(tax.table)
rownames(tax.table) <- data$Feature.ID
tax.table <- tax_table(tax.table)


# make into phyloseq object 
ps <- phyloseq(otu.table, tax.table, sam.data.ps)
ps <- subset_taxa(ps, Kingdom !="NA") # removes unknown kingdoms 
ps <- subset_taxa(ps, Kingdom !="") # removes blank kingdoms
ps <- subset_taxa(ps, Kingdom !="Unassigned") # removes unassigned kingdom 
ps <- subset_taxa(ps, Phylum != "") # removes unknown phyla 
ps.genus <- tax_glom(ps, "Genus") # may be used for downstream analyses, this groups the phyloseq file by Genus, we may not use this and it takes a while to run so ONLY RUN IF NECESSARY


# make into ampvis object to calculate alpha diversity, this is possible in phyloseq but using ampvis creates a datafram that is easier to use (imo)
d <- amp_load(otutable = ps@otu_table,
              metadata = metadata,
              taxonomy = ps@tax_table) # you will get a warning here "could not find a column..." ignore this, it is okay

# make into microeco object 
df <- phyloseq2meco(physeq = ps)
df
df$filter_pollution(taxa = c("mitochondria", "chloroplast"))
df$sample_sums() %>% range

set.seed(123456)

t1 <- trans_abund$new(dataset = df, taxrank = "Phylum", ntaxa = 10, groupmean = "Pipeline") # make plot with phylum, grouped by whatever variable 

t1$plot_bar() # use color_values to set your own colors, there are other adjustments you can use to make sure the bar plot is showing everthing you want to 

# use ampvis to make a heatmap 
amp_heatmap(data = d,
            tax_aggregate = "Genus", # show genera
                       tax_add = "Phylum", # show the phyla 
                       group_by = c("Treatment"), # group by categorical variables
                       tax_show = 25, # show 25 taxa 
                       plot_colorscale = "sqrt",
                       color_vector = c("#4A708B","#EEB422",  "#CD5B45"), # these can be changed
                       plot_values = FALSE, 
                       plot_legendbreaks = c(0,25,50),
                       facet_by = c("Pipeline")) + # facet by categorical variable
                      theme_bw() + 
                      theme(text = element_text(size =18),
                            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# use microeco to make venn diagrams 

dataset2 <- clone(df)
dataset2$tidy_dataset()
dataset2 <- dataset2$merge_samples(use_group = "Pipeline") # choose group
t1 <- trans_venn$new(dataset2, ratio = "numratio") # the numratio makes it so all the percentages add up to 100%, the default option is harder to interpret
t1$plot_venn() # use color_circle =  to set colors

