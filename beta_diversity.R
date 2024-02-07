# Beta Diversity 

# you should have already calculated alpha diversity, so the beginning of the script should look familiar 

# load libraries 

library(tidyverse)
library(phyloseq)
library(microeco)
library(file2meco)
library(vegan)
library(ampvis2)
library(ggpubr)
library(rstatix)

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

# make into microeco object 
df <- phyloseq2meco(physeq = ps)
df
df$filter_pollution(taxa = c("mitochondria", "chloroplast"))
df$sample_sums() %>% range
set.seed(123456)

df$cal_betadiv() # calculate beta diversity (diversity between groups)
t1 <- trans_beta$new(dataset = df, measure = "bray", group = "Location") # for group, choose the variable you are interseted in 
t1$cal_manova(manova_all = TRUE)
t1$res_manova # location is significant
t1 <- trans_beta$new(dataset = df, measure = "bray", group = "Plant") # you can run this multiple times to look at each categorical variable 
t1$cal_manova(manova_all = TRUE)
t1$res_manova # plant is not significant
t1$cal_anosim(group = "Treatment") # the anosim looks for differences *within* a group, if there is a lot of variablility within a group, this is important to consider when interpreting the data, in this case it is better if the p-value is > 0.05
t1$res_anosim


# this uses all ASVs but, we can also use the ps.genus that we made before and compare the two, sometimes, there are interesting patterns when we just use Genus because those are the seuqneces that we have more taxonomic information about 

# now we can plot it, we have to look at the data to determine if a NMDS or PCOA is better (I am guessing PCOA)
t1$cal_ordination(ordination = "PCoA")
t1$plot_ordination(plot_type = c("point"), 
                   plot_color = "Pipeline", # choose the variables for color and shape
                   plot_shape =  "Replicate", point_size = 5) + 
  theme(panel.grid = element_blank()) + 
  geom_vline(xintercept = 0, linetype = 2) + # adds a grid at 0,0, this makes it easier to see differences
  geom_hline(yintercept = 0, linetype = 2) + 
  theme_bw()
  # scale_color_manual(values = pal) + # optional change to color pallette (recommended)
  # facet_grid(.~Treatment) # optional facet by additional variable 

# on the microeco website, there are lots of variation for the POCA plot that we can look at, this is a pretty basic one 
