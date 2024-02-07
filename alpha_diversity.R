#  Alpha Diversity 

# this code will make plots and calculate statistics for alpha diversity (richness and diversity)

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


# make dataframe for alpha diversity 
alpha <- amp_alpha_diversity(d)

# create function to summarize data, this is used later to create error bars on the bar plots 

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
} # you only have to run this once 

# make dataframes for observed and shannon diversity (we can add others but these are nice ones to start with)

# you may be asked to install the package "plyr" 

observed.16S <- data_summary(alpha, varname="ObservedOTUs", 
                             groupnames=c("Pipeline", "Treatment"))
shannon.16S <- data_summary(alpha, varname="Shannon", 
                            groupnames=c("Pipeline", "Treatment"))

# make dataframe with the statistics, this will be used in the bar plot 
stat.test.observed <- alpha %>%
  group_by(Treatment) %>% #change to whatever variable you are interested in grouping by
  t_test(ObservedOTUs ~ Pipeline) %>% # change formula to match the variable you are interested in, t-tests compare two groups and works well for the bar plot we are about to create but an anova would work too, we will use and ANOVA later. The formula is (y ~ x) so we might do (ObservedOTUs ~ Day)
  adjust_pvalue(method = "bonferroni") %>% # since we have multiple samples, a p-adjustment is necessary 
  add_significance() # adds category to denote significance (P<0.05)

stat.test.observed <- stat.test.observed %>% add_xy_position(x = "Pipeline") # adding the xy postition allows us to plot it

observed.16s.plot <- ggplot(observed.16S, 
                            aes(x=Pipeline, y=ObservedOTUs)) + # we made "observered.16S above, x and y are the same as the formula in the stat.test
  geom_bar(stat="identity", color="black", 
           position=position_dodge(), 
           fill = "#CD5B45") + # change to whatever color you want 
  geom_errorbar(aes(ymin=ObservedOTUs-sd, # + - the standard deviation
                    ymax=ObservedOTUs+sd), width=.2, 
                position=position_dodge(.9)) +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size =12)) + 
  facet_grid(.~Treatment) + # facet by category (optional)
  stat_pvalue_manual(stat.test.observed) # add p-value 

# view plot

observed.16s.plot


# repeat with shannon 

stat.test.shannon <- alpha %>%
  group_by(Treatment) %>%
  t_test(Shannon ~ Pipeline) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test.shannon <- stat.test.shannon %>% add_xy_position(x = "Pipeline")

shannon.16S.plot <- ggplot(shannon.16S, aes(x=Pipeline, y=Shannon)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge(), fill = "#CD5B45") +
  geom_errorbar(aes(ymin=Shannon-sd, ymax=Shannon+sd), width=.2, position=position_dodge(.9)) + theme_bw() + theme(legend.position = "none", text = element_text(size =12)) + facet_grid(.~Treatment) + stat_pvalue_manual(stat.test.shannon)

shannon.16S.plot

# now we can test the difference statistsically, we already did a t-test but there are more accurate and comprehensive methods of testing the differences between communities

# make into microeco object 
df <- phyloseq2meco(physeq = ps)
df
df$filter_pollution(taxa = c("mitochondria", "chloroplast")) # remove mitochondria and chloroplast because aren't "true bacteria"
df$sample_sums() %>% range # look at the range of the data 

# before calculating anything, set the seed -- this tells R to use the same random numbers when cacluating the values -- if you don't set this then everytime you run the code you will get  *slightly* different values 

set.seed(12345)

t1 <- trans_alpha$new(dataset = df) # creates microeco trans_alpha object


## multi factor ANOVA ## 
t1$cal_diff(method = "anova", formula = c("Pipeline", "Location", "Type"), anova_post_test = "HSD.test") # multifactor ANVOA, here in the "Formula" you can enter all the variables you are interested in 
head(t1$res_diff)
View(t1$res_diff) # view the results in a new window, hopefully we will see some significant differences ;) 


# now you have successfully calculated and ploted alpha diveristy, there are other ways we can approach it too -- we can make box plots or facet it differently, it all depends on how the data looks. 
