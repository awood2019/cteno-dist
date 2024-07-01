#analyzing ctenophore GoMx eDNA sequences by using Allen's ctenophore sequences as a reference

#load packages
library(tidyverse)
library(dada2)
library(viridis) ; packageVersion("viridis") 

#set working directory
setwd("~/Desktop/AW RStudio/data/gomx-cteno-dist")

#--------LULU to curate ASVs---------

#install packages needed for lulu
install.packages("devtools")
install_github("tobiasgf/lulu")

#load packages
library(devtools)
library(lulu)

#read in files
gomx.table <- read.delim('table.tsv') #read in GoMx ASVs with their counts for each sample
gomx.seqids <- read.csv('rep-seqs-phylum.csv') #read in BLAST results of GoMx ASVs
matchlist <- read.delim('match_list.txt', header = FALSE) #read in match list of ASVs BLASTed against themselves

#put files in format needed for LULU by changing sequences to ASV1, ASV2, etc
colnames(gomx.table)[1] = "Sequence" #rename column of ASVs to "Sequence"
gomx.joined.table <- left_join(gomx.seqids, gomx.table, by = "Sequence") #merge sample count table with BLAST result table 
gomx.joined.table <- subset(gomx.joined.table, select = -c(Sequence, Phylum)) #remove columns we don't need for LULU 
gomx.final.table <- subset(gomx.joined.table, select = -c(ASV)) #remove column of ASV sequences since all columns must be sample counts
row.names(gomx.final.table) <- gomx.joined.table$ASV #make the rownames of the count table the ASV sequences to fit LULU requirements

#run LULU
curated.table <- lulu(gomx.final.table, matchlist, minimum_match = 95) #with 95% minimum match between sequences
curated.table.df <- curated.table$curated_table #assign table of curated ASVs to dataframe

#join curated ASV table with table of sequence IDs
curated.table.df$ASV <- row.names(curated.table.df) #assign row names of curated table to new column
curated.table.df <- left_join(curated.table.df, gomx.seqids, by = "ASV") #join curated df with sequence ID df by ASV column
row.names(curated.table.df) <- curated.table.df$ASV #make ASV ID the row names again since joining got rid of them

#select only ctenophores for further analysis
cteno.curated.table <- curated.table.df %>%
  subset(Phylum == "Ctenophora") 
  

#-------classify the curated ASVs--------

#copy sequence IDs to vector to compare to reference files
seqs <- cteno.curated.table$Sequence

#assignTaxonomy using FASTA file as reference
taxa <- assignTaxonomy(seqs, "Anth-28S-eDNA_ctenos_with_outgroups_assignTaxonomy.fasta", multi=TRUE, minBoot = 80,
                       taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) #classifying to the genus level with assignTaxonomy
unname(taxa)

#assignSpecies using FASTA file as reference
species <- addSpecies(taxa, "Anth-28S-eDNA_ctenos_with_outgroups_assignSpecies.fasta", allowMultiple=TRUE) #finding 100% matches to our reference database of Gulf of Mexico ctenos with assignSpecies
unname(species)
unique(species [, 7]) #see how many were identified to species level

#-------create summary table of LULU curated ASVs-------

#convert from vector to dataframe
cteno.species.df <- as.data.frame(species)

#import table of uncurated cteno ASVs so we can name rows of curated table to be the same
setwd("~/Desktop/AW RStudio/results/gomx-cteno-dist")
uncurated.classifier <- read_tsv("Cteno_classifier_results.tsv")

#join table of curated ASVs with uncurated ASVs to get the sequence labels
curated.and.uncurated <- left_join(uncurated.classifier, cteno.species.df, by = c("seq" = "Sequence"))
curated.and.uncurated$LULUresults <- ifelse(is.na(curated.and.uncurated$Name.y), "discarded", "curated") #create column for whether ASV was discarded by LULU

cteno.species.df <- curated.and.uncurated[,-2:-8] #remove classification columns from uncurated classifier
cteno.species.df <- subset(cteno.species.df, select = -c(totalcount, Name.y)) #remove redundant name column 

#save as a table
write.table(cteno.species.df, file = 'Cteno_classifier_lulu_results.tsv', sep = "\t", row.names = FALSE, quote=FALSE) #writing the ASV counts table with the taxonomic classifications of each cteno ASV


#-------DON'T NEED ANYMORE? phyloseq code adapted from GoMx phyla distribution script------

#load packages
library(tidyverse) ; packageVersion("tidyverse") 
library(phyloseq) ; packageVersion("phyloseq") 
library(vegan) ; packageVersion("vegan") 
library(DESeq2) ; packageVersion("DESeq2") 
library(dendextend) ; packageVersion("dendextend") 
library(viridis) ; packageVersion("viridis") 
library("ggplot2")

#create components of phyloseq object: taxonomy table, count table, and sample data table.
tax_tab <- cteno.species.df #taxonomy table with ASV sequences as row names and taxonomic classification
count_tab <- cteno.curated.table #count table with ASV counts for each sample
sample_info_tab <- read_csv("anth-28S-sampledata_20231016.csv") #load sample data table with sample metadata

#coerce tables into proper format to make phyloseq object
#tax_tab_phy: includes taxonomic information for each representative (ASV) sequence
tax_tab_phy <- as.matrix(tax_tab)

#count_tab_phy: includes all ASVs and their abundances in each sample (row.names must match row.names of tax_tab_phy)
row.names(count_tab) <- count_tab$Sequence
count_tab <- subset(count_tab, select = -c(ASV, Sequence, Phylum))
count_tab_phy <- as.matrix(count_tab)

#sample_info_tab_phy: table that includes sample information for all samples (row.names must equal col.names in count table)
sample_info_tab <- sample_info_tab %>% mutate(depth_bin = cut_width(sample_info_tab$Depth, width = 10, boundary = 0)) #create column for depth range as a factor
sample_info_tab_phy <- sample_info_tab
sample_info_tab_phy <- sample_info_tab_phy[-c(55,56),] #delete the last 2 rows because they have NAs across the board
sample_data <- sample_data(sample_info_tab_phy) #convert to phyloseq component now because row names get changed by sample_data command
row.names(sample_data) <- sample_data$File.name #change row names to match file name

#make phyloseq object with just the count table and taxonomy table
ASV_physeq <- phyloseq(otu_table(count_tab_phy, taxa_are_rows = TRUE), tax_table(tax_tab_phy), sample_data)
ASV_physeq <- prune_taxa(taxa_sums(ASV_physeq) > 0, ASV_physeq) #pruning out ASVs with zero counts
saveRDS(ASV_physeq, 'allphy_physeq.rds') #save phyloseq object

#transform phyloseq object to dataframe for easy viewing
df_ASV_physeq <- ASV_physeq %>% psmelt() #melt phyloseq object to long dataframe
head(df_ASV_physeq)

