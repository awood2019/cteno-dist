#analyzing ctenophore GoMx eDNA sequences by using Allen's ctenophore sequences as a reference

#load packages
BiocManager::install("dada2")
library(dada2)

#--------code from GoMx phyla distribution script to make phyloseq object combining taxonomy, sample metadata, and abundance data-------

#load packages
library(tidyverse) ; packageVersion("tidyverse") 
library(phyloseq) ; packageVersion("phyloseq") 
library(vegan) ; packageVersion("vegan") 
library(DESeq2) ; packageVersion("DESeq2") 
library(dendextend) ; packageVersion("dendextend") 
library(viridis) ; packageVersion("viridis") 
library("ggplot2")

#set working directory 
setwd("~/Desktop/AW RStudio/data/gomx-phy-dist")

#load components of phyloseq object: taxonomy table, count table, and sample data table.
tax_tab <- read_csv("rep-seqs-phylum.csv", show_col_types = FALSE) #loading taxonomy table w/ ASVs, sequence, & phyla
count_tab <- read_delim("table.tsv") #loading count table w/ ASV counts for each sample
sample_info_tab <- read_csv("anth-28S-sampledata_20231016.csv") #loading sample data table w/ sample metadata

#coerce tables into proper format to make phyloseq object
#tax_tab_phy: includes taxonomic information for each representative (ASV) sequence
phylum <- tax_tab$Phylum #pulling out phylum column from taxonomy table
tax_tab_phy <- tibble(phylum) #making phyla into a tibble containing phylum for each sequence
tax_tab_phy <- as.matrix(tax_tab_phy) #make tibble into matrix
row.names(tax_tab_phy) <- tax_tab$Sequence #make sequence column the row names
tax_tab_phy[is.na(tax_tab_phy)] <- "< 85% similarity to top BLAST hit" #change NA values to more accurate description

#count_tab_phy: includes all ASVs and their abundances in each sample (row.names must match row.names of tax_tab_phy)
count_tab_phy <- select(count_tab, -"...1") #delete this weird column
row.names(count_tab_phy) <- count_tab$...1 #make sequences the row names (ignore warning message)

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
setwd("~/Desktop/AW RStudio/results/gomx-phy-dist")
write.table(df_ASV_physeq, file = "GoM_allASVs.tsv", sep = "\t",row.names = FALSE, quote = FALSE) #save file as a .tsv

#----------end of the code from GoM phyla distribution script--------

#-------------assign taxonomy to ctenophore ASVs--------------

#switch working directory
setwd("~/Desktop/AW RStudio/data/gomx-cteno-dist")

#select only the ctenophores
cteno.table <- df_ASV_physeq %>%
  subset(phylum == "Ctenophora") %>%
  filter(Abundance > 0) #only the samples that contained ctenophores
write.table(cteno.table, file = "Cteno_preclassifier.tsv", sep = "\t",row.names = FALSE, quote = FALSE) #save file as a csv

#create vector of only unique sequences from cteno.table so we can compare them to the FASTA files
seqs <- unique(cteno.table$OTU)

#assignTaxonomy with Allen's new FASTA file 
taxa <- assignTaxonomy(seqs, "AC_Cteno-28S-assignTaxonomy.fasta", multi=TRUE, minBoot = 80,
                       taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) #classifying to the genus level with assignTaxonomy
unname(taxa)

#assignSpecies using Allen's new FASTA file 
taxa.species <- addSpecies(taxa, "AC_Cteno-28S-assignSpecies.fasta", allowMultiple=TRUE) #finding 100% matches to our reference database of Gulf of Mexico ctenos with assignSpecies
unique(taxa.species [, 7]) #see how many were identified to species level

#convert from vector to dataframe
cteno.species.df <- as.data.frame(taxa.species)
cteno.species.df$seq <- row.names(cteno.species.df)

#add abundance of each unique ASV to dataframe
length(unique(cteno.species.df[["seq"]])) #find how many unique ASVs are in table

counts_per_ASV <- cteno.table %>% #make table of abundance counts for each unique ctenophore ASV
  group_by(OTU) %>%
  summarize(totalcount = sum(Abundance)) 
names(counts_per_ASV)[names(counts_per_ASV) == "OTU"] <- "seq" #name this column seq so we can combine it with cteno.species.df
names(counts_per_ASV)[names(counts_per_ASV) == "totalcount"] <- "unique abundance" #name this column so we know it represents the counts of unique ASVs present in dataset

#combine abundance counts table of unique ctenophores w/ classifier table
cteno.species.df <- left_join(cteno.species.df, counts_per_ASV, by = "seq")

#save as a table
setwd("~/Desktop/AW RStudio/results/gomx-cteno-dist")
cteno.species.df[is.na(cteno.species.df$'Family')] <- "not identified to Family level" #NEED TO FIX change NA values to more accurate description
write.table(cteno.species.df, file = 'Cteno_classifier_results.tsv', sep = "\t", row.names = FALSE, quote=FALSE) #writing the ASV counts table with the taxonomic classifications of each cteno ASV

#create table from cteno.table with only depth, method & abundance 
cteno.depth <- cteno.table %>%
  select(OTU, Abundance, depth_bin, CTD.ROV)
cteno.species.depth <- right_join(cteno.species.df, cteno.depth, by = c("seq" = "OTU")) #join with table of cteno species taxonomy
cteno.species.depth <- cteno.species.depth %>% #take out all samples that were blanks or controls
  filter(CTD.ROV != "Negative" & CTD.ROV != "NTC")

#figure out which depth ranges were sampled at so we can make depth bins for plots
unique(cnid.species.depth$depth_bin)

#plot by genus
ggplot(cteno.species.depth, aes(x=factor(depth_bin, level=c('[0,10]', '(40,50]', '(50,60]', '(60,70]', '(70,80]', '(80,90]', '(110,120]', '(440,450]', '(450,460]','(460,470]', '(470,480]', '(520,530]', '(530,540]' )), y = Abundance, fill = Genus)) + #x-axis = depth, y-axis = ASV abundance - plotted by genus
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(.~CTD.ROV, scale = "free_x", space = "free_x") +
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 

#plot by family
ggplot(cteno.species.depth, aes(x=factor(depth_bin, level=c('[0,10]', '(40,50]', '(50,60]', '(60,70]', '(70,80]', '(80,90]', '(110,120]', '(440,450]', '(450,460]','(460,470]', '(470,480]', '(520,530]', '(530,540]' )), y = Abundance, fill = Family)) + #x-axis = depth, y-axis = ASV abundance - plotted by genus
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(.~CTD.ROV, scale = "free_x", space = "free_x") + #facet by sampling method
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  ggtitle("Ctenophora Family Distribution in Gulf of Mexico")+
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 

#plot by order
ggplot(cteno.species.depth, aes(x=factor(depth_bin, level=c('[0,10]', '(40,50]', '(50,60]', '(60,70]', '(70,80]', '(80,90]', '(110,120]', '(440,450]', '(450,460]','(460,470]', '(470,480]', '(520,530]', '(530,540]' )), y = Abundance, fill = Order)) + #x-axis = depth, y-axis = ASV abundance - plotted by genus
  geom_bar(position="fill", stat = "identity") + #position=fill graphs abundance as a proportion out of the total, stat=identity tells ggplot to calculate sum of the y var grouped by the x var
  facet_grid(.~CTD.ROV, scale = "free_x", space = "free_x") +
  ylab("Percentage of ASVs recovered") + 
  xlab("Depth (m)") +
  theme_classic() +
  scale_fill_viridis(discrete=TRUE, option="turbo") 

na.order <- subset(cteno.species.depth, is.na(cteno.species.depth$Order))
na.order[1,8]
na.order[2,8]

#---------------plotting ctenophore ASVs in PCoA plot with UniFrac distance matrix---------------

#load packages needed for PCoA ordination plot
library(vegan)
library(dplyr)
library(ape)
library(phyloseq)
library(Biostrings)
library(ShortRead)
library(tidysq)
library(ggalt)
library(phytools)

#set working directory
setwd("~/Desktop/AW RStudio/data/gomx-cteno-dist")

#make phyloseq object with only ctenophores
cteno_physeq = subset_taxa(ASV_physeq, phylum == "Ctenophora") 

#take out all negative and blank samples
cteno_physeq_pos <- subset_samples(cteno_physeq, CTD.ROV != "Negative" & CTD.ROV !="NTC")

#remove rows with 0 ctenophores so the distance matrix can be calculated
cteno_physeq_0 <- prune_samples(sample_sums(cteno_physeq_pos) >0, cteno_physeq_pos) 

#also remove rows with less than 10 and 100 ctenophores to see how ordination changes 
cteno_physeq_10 <- prune_samples(sample_sums(cteno_physeq) >10, cteno_physeq) 
cteno_physeq_100 <- prune_samples(sample_sums(cteno_physeq) >100, cteno_physeq) 

#export table of ASVs so we can put it into Geneious to make a phylogenetic tree
tax_cteno <- tax_table(cteno_physeq_0) #take component of ASVs out of phyloseq object
write.csv(tax_cteno, "~/Desktop/AW RStudio/tax_cteno.csv") #save as csv file

#add tree of ctenophores to each phyloseq object
cteno_tree <- read.nexus("/Users/quattrinia/Desktop/AW RStudio/data/gomx-cteno-dist/tax_cteno new alignment FastTree Tree.nex") #import Geneious tree of aligned ctenophores
rooted_tree <- midpoint.root(cteno_tree)
physeq_tree_abdndce0 <- merge_phyloseq(cteno_physeq_0, rooted_tree) #merge tree with existing phyloseq object for abundance >0
physeq_tree_abdndce10 <- merge_phyloseq(cteno_physeq_10, rooted_tree) #merge tree with existing phyloseq object for abundance >10
physeq_tree_abdndce100 <- merge_phyloseq(cteno_physeq_100, rooted_tree) #merge tree with existing phyloseq object for abundance >100

#calculate Unifrac distance matrix to each phyloseq object
uni_matrix_0 <- UniFrac(physeq_tree_abdndce0, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE) #distance matrix for phyloseq object with abundance >0
uni_matrix_10 <- UniFrac(physeq_tree_abdndce10, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE) #distance matrix for phyloseq object with abundance >10
uni_matrix_100 <- UniFrac(physeq_tree_abdndce100, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE) #distance matrix for phyloseq object with abundance >100

#perform ordination for each phyloseq object
cteno_ordination_0 <- ordinate(physeq_tree_abdndce0, method = "PCoA", uni_matrix_0, weighted=TRUE) #ordination for phyloseq object with abundance >0
cteno_ordination_10 <- ordinate(physeq_tree_abdndce10, method = "PCoA", uni_matrix_10, weighted=TRUE) #ordination for phyloseq object with abundance >10
cteno_ordination_100 <- ordinate(physeq_tree_abdndce100, method = "PCoA", uni_matrix_100, weighted=TRUE) #ordination for phyloseq object with abundance >100

#make PCOA plot for each phyloseq object
#plot for phyloseq object with abundance >0
plot_ordination(physeq_tree_abdndce0, cteno_ordination_0, color = "Depth", shape="CTD.ROV") + #define the point color by depth and point shape by sampling method 
  geom_encircle(aes(group=Site, fill=Site), alpha=0.2, s_shape=1, expand=0) + #add polygon shapes to enclose points by site
  theme_classic()+
  scale_color_continuous(type="viridis", option="D", direction=-1) +
  ggtitle("PCoA of Ctenophore ASVs with Abundance >0")

#plot for phyloseq object with abundance >10
plot_ordination(physeq_tree_abdndce10, cteno_ordination_10, color = "Depth", shape="CTD.ROV") + #define the point color by depth and point shape by sampling method 
  geom_encircle(aes(group=Site, fill=Site), alpha=0.2, s_shape=1, expand=0) + #add polygon shapes to enclose points by site
  theme_classic()+
  scale_color_continuous(type="viridis", option="D", direction=-1) +
  ggtitle("PCoA of Ctenophore ASVs with Abundance >10")

#plot for phyloseq object with abundance >100
plot_ordination(physeq_tree_abdndce100, cteno_ordination_100, color = "Depth", shape="CTD.ROV") + #define the point color by depth and point shape by sampling method 
  geom_encircle(aes(group=Site, fill=Site), alpha=0.2, s_shape=1, expand=0) + #add polygon shapes to enclose points by site
  theme_classic()+
  scale_color_continuous(type="viridis", option="D", direction=-1) +
  ggtitle("PCoA of Ctenophore ASVs with Abundance >100")

#trying Cori's code for % of variance explained
Axis1 <- cteno_ordination_0$vectors[,1]
Axis2 <- cteno_ordination_0$vectors[,2]
Var1 <- round(cteno_ordination_0$values[1,2], digits=4)*100
Var2 <- round(cteno_ordination_0$values[2,2], digits=4)*100
X_axis <- paste("Axis1"," ","[",Var1,"%]",sep="")
Y_axis <- paste("Axis2"," ","[",Var2,"%]",sep="")

plot_ordination(physeq_tree_abdndce0, cteno_ordination_0, color = "Depth", shape="CTD.ROV", label = X_axis) + #define the point color by depth and point shape by sampling method 
  geom_encircle(aes(group=Site, fill=Site), alpha=0.2, s_shape=1, expand=0) + #add polygon shapes to enclose points by site
  theme_classic()+
  scale_color_continuous(type="viridis", option="D", direction=-1) +
  labs(x=X_axis, y=Y_axis) +
  ggtitle("PCoA of Ctenophore ASVs with Abundance >0")

#---------------plotting ctenophore ASVs in PCoA plot with Bray-Curtis distance matrix---------------
cteno.table <- df_ASV_physeq %>%
  subset(phylum == "Ctenophora") %>% #only the samples that contained ctenophores
  subset(CTD.ROV != "Negative" & CTD.ROV !="NTC") %>%
  filter(Abundance > 5) 

cteno.abundance <- cteno.table %>% #make a table from cteno.table
  select(c(1:3)) %>% #select only the first 3 columns
  spread(Sample, Abundance, fill = 0) #group the table so that Sample name is column names and Abundance is the cell values

rownames(cteno.abundance) <- cteno.abundance$OTU #make the row names the OTU names

cteno.matrix <- t(as.matrix(select(cteno.abundance, -1))) #make a matrix with the OTU names as column names and samples as row names

sapply(cteno.matrix, class) #see which class each cell in the matrix is because we need them to all be numeric for distance matrix function

colSums(cteno.matrix) #see column sums for each matrix column to make sure none of them are 0
rowSums(cteno.matrix) #see row sums for each matrix row
any(is.na(cteno.matrix)) #check if there are any NA values in matrix

#calculate a distance matrix
dist <- vegdist(cteno.matrix, method = "bray") #use Bray-Curtis distance metric

rel <- decostand(dist, "total", 1)

#conduct PCoA (principal coordinate analysis)
PCOA <- pcoa(dist)
biplot.pcoa(PCOA) #plot PCoA

#make table of PCoA points and metadata to use in ggplot
PCOAinfo <- data.frame(PCOA$vectors[,1], PCOA$vectors[,2]) #pull out information from PCoA object that we want for ggplot
PCOAinfo$File.name <- rownames(PCOAinfo) #make column for file name so we can join with metadata table
metadata <- read_csv("anth-28S-sampledata_20231016.csv") #loading sample data table w/ sample metadata
ggplot.table <- left_join(metadata, PCOAinfo, by = 'File.name') #combine tables
ggplot.table <- ggplot.table %>% subset(`CTD/ROV` != "Negative" & `CTD/ROV` !="NTC") #remove negatives and blanks
colnames(ggplot.table) #get column names for table of ggplot info
ggplot.table$pcoa1 <- ggplot.table$PCOA.vectors...1. #rename axis 1 column 
ggplot.table$pcoa2 <- ggplot.table$PCOA.vectors...2. #rename axis 2 column
ggplot.table <- select(ggplot.table, -17, -18)

#plot PCoA using ggplot
ggplot(ggplot.table, aes(pcoa1, pcoa2, color = Depth, shape = `CTD/ROV`, label = `Sample name`)) +
  geom_point() +
  theme_classic() +
  scale_color_viridis_c(direction = -1) +
  geom_label() +
  geom_text(check_overlap = TRUE)

#trying Cori's code to see how % of variance explained by PCoA
axis1 <- PCOA$vectors[,1]
axis2 <- PCOA$vectors[,2]
var1 <- round(PCOA$values[1,2], digits=4)*100
var2 <- round(PCOA$values[2,2], digits=4)*100
x_axis <- paste("axis1"," ","[",var1,"%]", sep="")
y_axis <- paste("axis2"," ","[",var2,"%]",sep="")


p <- ggplot(ggplot.table, aes(axis1, axis2, color=Depth, shape=`CTD/ROV`, label = `Sample name`)) + 
  geom_point()
  theme_classic() + 
  labs(x = x_axis, y = y_axis)
p



