###Differential Abundance with DESeq2

Adapted from https://joey711.github.io/phyloseq-extensions/DESeq2.html
setwd("~/Desktop/eunice/Thesis/Qiime/Samples/OnlySamples/Filtered/DESeq/")

rm(list = ls ())


library("DESeq2")
library(dplyr)
library(tidyr)
library(ape)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(plotly)
library(tidyr)
library(naniar)
library(zoo)
library(lubridate)
library(qiime2R)

#OTU table (shared file)
#The OTU table as exported from qiime has a pound sign before the header row. You need to delete that pound sign in a text editor.
metadata <- read.csv("DESeqmetadata.csv", na.strings = c("","NA"), header=TRUE)
#metadata <- metadata2[-1,]
str(metadata)
metadata$BRD <- factor(metadata$BRD) 
order_groups <- metadata$ID
row.names(metadata) = metadata[,1]
metadata = metadata[,-1]


ASVs <- read_qza("~/Desktop/eunice/Thesis/Qiime/Samples/OnlySamples/Filtered/Qiime/table-filtered2.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #18010 ASVs
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
ASV_table <- t(ASV_table)

#Taxonomy of each OTU
##Adding taxonomy
#Taxonomy of each OTU
tax <- read_qza("~/Desktop/eunice/Thesis/Qiime/Samples/OnlySamples/Filtered/Qiime/taxonomy.qza")
tax <- as.data.frame(tax$data)
tax2 = separate(tax, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#This warning means that some cells are empty and that R is replacing empty cells with NA. Now there are other cells that are unclassified that  say, for example `s__` or `g__`. All these need to be removed and replaced with `NA`. 
#All this is OK except that in future use of the taxonomy table, these ASVs will be ignored because they are not classified. Why are ASVs not classified? Its because there is not a close enough match in the database. Just because there is not a good match in the database does not mean they don’t exist, so I wanted to make sure this data was not lost. So in my new code, from lines 300 – 334 I make it so that ASVs that are unclassified at any level are classified as the lowest taxonomic level for which there is a classification.

#All the strings that need to be removed and replaced with NA
na_strings <- c(" s__", " g__", " f__", " o__", " c__")

tax3 = replace_with_na_all(tax2, condition = ~.x %in% na_strings)

#This code is great because an ASV that is unclassified at a certain level are all listed as `NA`.
#Unfortunately this command changed ou Feature.ID names

#Next, all these `NA` classifications with the last level that was classified
tax3[] <- t(apply(tax3, 1, zoo::na.locf))
tax3 <- as.data.frame(tax3)
row.names(tax3) <- tax3[,1]
tax3 = tax3[,-c(1:2)]
tax.clean <- as.data.frame(tax3)
tax.clean$OTUs <- rownames(tax.clean)
#Would be good to check here to make sure the order of the two data frames was the same. You should do this on your own.

###Remove all the OTUs that don't occur in our OTU.clean data set
tax.final = tax.clean[row.names(tax.clean) %in% row.names(ASV_s),]

##Remove unneccesary information from the taxonomy names
tax.final$Phylum <- sub("D_0__*", "", tax.final[,1])
tax.final$Phylum <- sub("D_1__*", "", tax.final[,1])
tax.final$Class <- sub("D_0__*", "", tax.final[,2])
tax.final$Class <- sub("D_1__*", "", tax.final[,2])
tax.final$Class <- sub("D_2__*", "", tax.final[,2])
tax.final$Order <- sub("D_0__*", "", tax.final[,3])
tax.final$Order <- sub("D_1__*", "", tax.final[,3])
tax.final$Order <- sub("D_2__*", "", tax.final[,3])
tax.final$Order <- sub("D_3__*", "", tax.final[,3])
tax.final$Family <- sub("D_0__*", "", tax.final[,4])
tax.final$Family <- sub("D_1__*", "", tax.final[,4])
tax.final$Family <- sub("D_2__*", "", tax.final[,4])
tax.final$Family <- sub("D_3__*", "", tax.final[,4])
tax.final$Family <- sub("D_4__*", "", tax.final[,4])
tax.final$Genus <- sub("D_0__*", "", tax.final[,5])
tax.final$Genus <- sub("D_1__*", "", tax.final[,5])
tax.final$Genus <- sub("D_2__*", "", tax.final[,5])
tax.final$Genus <- sub("D_3__*", "", tax.final[,5])
tax.final$Genus <- sub("D_4__*", "", tax.final[,5])
tax.final$Genus <- sub("D_5__*", "", tax.final[,5])
tax.final$Species <- sub("D_0__*", "", tax.final[,6])
tax.final$Species <- sub("D_1__*", "", tax.final[,6])
tax.final$Species <- sub("D_2__*", "", tax.final[,6])
tax.final$Species <- sub("D_3__*", "", tax.final[,6])
tax.final$Species <- sub("D_4__*", "", tax.final[,6])
tax.final$Species <- sub("D_5__*", "", tax.final[,6])
tax.final$Species <- sub("D_6__*", "", tax.final[,6])

TaxASV <- merge(tax.final, ASVkey, by.x = 0, by.y = "ASVstring")
row.names(TaxASV) <- TaxASV[,10]
TaxASV = TaxASV[,-c(1,10)]

### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq

colnames(tax_table(physeq_deseq))
## Filter any non-baxteria, chloroplast and mitochondria
physeq_deseq %>%
  subset_taxa(Family != "Mitochondria" & 
                Genus != "Mitochondria" &
                Species != "Mitochondria" &
                Order != "Chloroplast" &
                Family != "Chloroplast" &
                Genus != "Chloroplast" &
                Species != "Chloroplast") -> physeq_deseq
physeq_deseq

#You need to run the phyloseq_to_df function
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]

NewASVTable <- prunetable[,c(1,10:140)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]
NewASVTable = t(NewASVTable)

### Checking how the data looks
## Make a plot to see the community in the two groups
# this prunes the taxa with abundance <2%
### CALCULATION OF THE ABUNDANCE OF EACH OTU  
otu.summary <- prop.table(as.matrix(NewASVTable), 1) 
str(otu.summary)
otu_abund <- colSums(otu.summary)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:17931)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)

#merging the abundance of each OTU with the metadata and the taxonomy file
str(metadata)
str(melt_otu)
meta_otu <- merge(metadata, melt_otu, by.x = 0, by.y = "Sample")
str(meta_otu)
meta_otu_tax <- merge(meta_otu, NewTax, by.x = "ASV", by.y = 0)
str(meta_otu_tax)
levels(meta_otu_tax$BRD)
str(metadata)
meta_otu_tax$Row.names <- factor(meta_otu_tax$Row.names, levels = order_groups)
summary(meta_otu_tax$Row.names) ###to check that all the samples have the same number of OTUs (6199 total, same value from the taxonomy file) 
meta_otu_tax$Family <- factor(meta_otu_tax$Family)
meta_otu_tax$Status <- factor(meta_otu_tax$BRD)
levels(meta_otu_tax$Status) <- list("Healthy"="0", "BRD"="1")
meta_otu_tax$Genus <- factor(meta_otu_tax$Genus)
meta_otu_tax$Phylum <- factor(meta_otu_tax$Phylum)
meta_otu_tax$ASV <- factor(meta_otu_tax$ASV)
str(meta_otu_tax)

### Checking the abundance of the most common taxa
## The abundance at a family level
family <- meta_otu_tax %>% 
  group_by(Family) %>% 
  summarise(Abundance = sum(Abundance))
attach(family)
family <- family[order(-Abundance),]

## Familly abundance only taking in consideration treatment effect in both days
family.2 <- meta_otu_tax %>% 
  group_by(Family, Status) %>% 
  summarise(Abundance = sum(Abundance))
attach(family.2)
family.2 <- family.2[order(-Abundance),]

### Abundance at a phlyum level
phylum <- meta_otu_tax %>% 
  group_by(Phylum) %>% 
  summarise(Abundance = sum(Abundance))
attach(phylum)
phylum <- phylum[order(-Abundance),]

## Phylum abundance only taking in consideration treatment effect in both days
phylum.2 <- meta_otu_tax %>% 
  group_by(Phylum, Status) %>% 
  summarise(Abundance = sum(Abundance))
attach(phylum.2)
phylum.2 <- phylum.2[order(-Abundance),]

## Phylum abundance only taking in consideration season effect in both days
phylum13.3 <- meta_otu_13 %>% 
  group_by(phylum, season) %>% 
  summarise(Abundance = sum(Abundance))
attach(phylum13.3)
phylum13.3 <- phylum13.3[order(-Abundance),]

### Abundance at a genus level
genus <- meta_otu_tax %>% 
  group_by(Genus) %>% 
  summarise(Abundance = sum(Abundance))
attach(genus)
genus <- genus[order(-Abundance),]

## Genus abundance only taking in consideration treatment effect
genus.2 <- meta_otu_tax %>% 
  group_by(Genus, Status) %>% 
  summarise(Abundance = sum(Abundance))
attach(genus.2)
genus.2 <- genus.2[order(-Abundance),]

## PHYLUM LEVEL
num_genera <- 97 # we need 100 OTUs in order to get the 25 most abundant Genus

melt_otu1 <- reshape2::melt(otu.summary_sorted[, c(1:num_genera)])
colnames(melt_otu1) <- c("Sample", "OTU", "Abundance")
tail(melt_otu1)

#Putting it all together: merge melt_otu, metadata, taxonomy tables
meta_otu1 <- merge(metadata, melt_otu1, by.x = 0, by.y = "Sample")
str(meta_otu1)
meta_otu_tax1 <- merge(meta_otu1, NewTax, by.x = "OTU", by.y = 0)
str(meta_otu_tax1)
meta_otu_tax1$Row.names <- factor(meta_otu_tax1$Row.names, levels = order_groups)
summary(meta_otu_tax1$Row.names) ###to check that all the samples have the same number of OTUs (346 total) 
meta_otu_tax1$Family <- factor(meta_otu_tax1$Family)
meta_otu_tax1$Status <- factor(meta_otu_tax1$BRD)
levels(meta_otu_tax1$Status) <- list("Healthy"="0", "BRD"="1")
meta_otu_tax1$Phylum <- factor(meta_otu_tax1$Phylum)
meta_otu_tax1$Genus <- factor(meta_otu_tax1$Genus)
meta_otu_tax1$Family <- factor(meta_otu_tax1$Family)
levels(meta_otu_tax1$Phylum)

##Whole phylum abuundance 
str(meta_otu_tax1)

### Calculation of the Phylum relative abundance for each time and treatment
PhylumAB <- meta_otu_tax1 %>% 
  group_by(Row.names, Status, Date.Collection, Phylum) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Status, Date.Collection, Phylum) %>%
  summarise(taxa.average = mean(taxa.sum)) ### relative abundance
str(PhylumAB)
PhylumAB$Phylum <- factor(PhylumAB$Phylum)
PhylumAB$Date.Collection <- factor(PhylumAB$Date.Collection)
levels(PhylumAB$Phylum)
levels(PhylumAB$Date.Collection)
levels(PhylumAB$Date.Collection) <- list("7/14"="7/14/20","7/21"="7/21/20","8/12"="8/12/20", "8/19"="8/19/20",
                                     "8/26"="8/26/20", "9/10"="9/10/20", "9/30"="9/30/20", "10/14"="10/14/20", 
                                     "10/28"="10/28/20","11/4"="11/4/20", "11/11"="11/11/20", "11/18"="11/18/20", "12/2"="12/2/20"
)


my_colors <- c(
  '#a6cee3','#1f78b4','#b3df8a','#33a03c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab3d6','#6a3d9a','#ffff99','#b15938', 
  "#CBD588", "#5F7FC7", "orange","#DA5734", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14385", "#653936", "#C84348", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

#Plot the graph 
ggplot(PhylumAB, aes(x = Date.Collection, y = taxa.average, fill =Phylum)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  facet_grid(Status~.)+
  #facet_wrap(vars(BRD), scales = "free") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=13)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 13), axis.text.y = element_text(color = "black", size = 13)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Relative Abundance Phylum (Top 8)")) +  labs(x='Date of Collection')

### FAMILY LEVEL
num_genera <- 44 # s

melt_otu3 <- reshape2::melt(otu.summary_sorted[, c(1:num_genera)])
colnames(melt_otu3) <- c("Sample", "OTU", "Abundance")
tail(melt_otu3)

#Putting it all together: merge melt_otu, metadata, taxonomy tables
meta_otu3 <- merge(metadata, melt_otu3, by.x = 0, by.y = "Sample")
meta_otu_tax3 <- merge(meta_otu3, NewTax, by.x = "OTU", by.y = 0)
str(meta_otu_tax3)
meta_otu_tax3$Row.names <- factor(meta_otu_tax3$Row.names, levels = order_groups)
summary(meta_otu_tax3$Row.names) ###to check that all the samples have the same number of OTUs (120) 
meta_otu_tax3$Family <- factor(meta_otu_tax3$Family)
meta_otu_tax3$Status <- factor(meta_otu_tax3$BRD)
levels(meta_otu_tax3$Status) <- list("Healthy"="0", "BRD"="1")
meta_otu_tax3$Phylum <- factor(meta_otu_tax3$Phylum)
meta_otu_tax3$Genus <- factor(meta_otu_tax3$Genus)
meta_otu_tax3$Family <- factor(meta_otu_tax3$Family)
levels(meta_otu_tax3$Phylum)

FamilyAB <- meta_otu_tax3 %>% 
  group_by(Row.names, Status, Date.Collection, Family) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Status, Date.Collection, , Family) %>%
  summarise(taxa.average = mean(taxa.sum)) 
str(FamilyAB)
FamilyAB$Family <- factor(FamilyAB$Family)
FamilyAB$Date.Collection <- factor(FamilyAB$Date.Collection)
levels(FamilyAB$Family)
levels(FamilyAB$Date.Collection)
levels(FamilyAB$Date.Collection) <- list("7/14"="7/14/20","7/21"="7/21/20","8/12"="8/12/20", "8/19"="8/19/20",
                                         "8/26"="8/26/20", "9/10"="9/10/20", "9/30"="9/30/20", "10/14"="10/14/20", 
                                         "10/28"="10/28/20","11/4"="11/4/20", "11/11"="11/11/20", "11/18"="11/18/20", "12/2"="12/2/20"
)

#Plot the graph 
ggplot(FamilyAB, aes(x = Date.Collection, y = taxa.average, fill =Family)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  facet_grid(Status~.)+
  ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.0, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  #theme(axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=13)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 13), axis.text.y = element_text(color = "black", size = 13)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Average Relative Abundance Family (Top 20)")) +  labs(x='Date of Collection')

#Genus level
num_genera <- 10 # we need 100 OTUs in order to get the 25 most abundant Genus

melt_otu4 <- reshape2::melt(otu.summary_sorted[, c(1:num_genera)])
colnames(melt_otu4) <- c("Sample", "OTU", "Abundance")
tail(melt_otu4)


#Putting it all together: merge melt_otu, metadata, taxonomy tables
meta_otu4 <- merge(metadata, melt_otu4, by.x = 0, by.y = "Sample")
meta_otu_tax4 <- merge(meta_otu4, NewTax, by.x = "OTU", by.y = 0)
str(meta_otu_tax4)
meta_otu_tax4$Row.names <- factor(meta_otu_tax4$Row.names, levels = order_groups)
summary(meta_otu_tax4$Row.names) ###to check that all the samples have the same number of OTUs (53) t)
meta_otu_tax4$Family <- factor(meta_otu_tax4$Family)
meta_otu_tax4$Status <- factor(meta_otu_tax4$BRD)
levels(meta_otu_tax4$Status) <- list("Healthy"="0", "BRD"="1")
meta_otu_tax4$Phylum <- factor(meta_otu_tax4$Phylum)
meta_otu_tax4$Genus <- factor(meta_otu_tax4$Genus)
meta_otu_tax4$Family <- factor(meta_otu_tax4$Family)
levels(meta_otu_tax3$Phylum)


##Whole family abuundance 
GenusAB <- meta_otu_tax4 %>% 
  group_by(Row.names, Status, Date.Collection, Genus) %>% 
  summarise(taxa.sum = sum(Abundance)) %>%
  group_by(Status, Date.Collection, Genus) %>%
  summarise(taxa.average = mean(taxa.sum)) 
str(PhylumAB)
GenusAB$Genus <- factor(GenusAB$Genus)
GenusAB$Date.Collection <- factor(GenusAB$Date.Collection)
levels(GenusAB$Genus)
levels(GenusAB$Date.Collection)
levels(GenusAB$Date.Collection) <- list("7/14"="7/14/20","7/21"="7/21/20","8/12"="8/12/20", "8/19"="8/19/20",
                                         "8/26"="8/26/20", "9/10"="9/10/20", "9/30"="9/30/20", "10/14"="10/14/20", 
                                         "10/28"="10/28/20","11/4"="11/4/20", "11/11"="11/11/20", "11/18"="11/18/20", "12/2"="12/2/20"
)

# PLOT FOR THE FIRST 25 GENUS
ggplot(GenusAB, aes(x = Date.Collection, y = taxa.average, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  scale_fill_manual(values = my_colors) +
  facet_grid(Status~.)+
  ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  #theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  #theme(axis.text.x = element_text(hjust = 1, vjust = 0.5)) +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .6, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  theme(strip.text = element_text(size = 13, face = "bold")) +
  theme(legend.text = element_text(size=13)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=13, face="bold"), axis.title.y = element_text(color="black", size=13, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 13), axis.text.y = element_text(color = "black", size = 13)) +
  #theme(text=element_text(family="Times New Roman")) +
  ylab(paste0("Average Relative Abundance Genus (Top 20)")) +  labs(x='Date of Collection')

## Running DESeq
## creating a new phyloseq object
#### DESEq differentially```
#To use DESeq, we need no zeros in our OTU table. So we will edit the table + 1

NewASVTable2 <- NewASVTable + 1

### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(NewASVTable2), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq

levels(metadata$BRD)

# establishing the model
diagdds = phyloseq_to_deseq2(physeq_deseq, ~ BRD)
#PenCode needs to be factor
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
head(diagdds)
resultsNames(diagdds)

my_contrast = c("BRD", "0", "1")
res = results(diagdds, contrast = my_contrast, cooksCutoff = FALSE, alpha=0.05)
summary(res)
res
res <- as.data.frame(res)
alpha = 0.05
#sigtab = res ### No significant results
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_deseq)[rownames(sigtab), ], "matrix"))
sigtab

sigtab$High_low <- ifelse(
  sigtab$log2FoldChange < -1.00, 'High in Healthy',
  ifelse(sigtab$log2FoldChange > 1.00, 'High in BRD',
         'Mid Change'))
#write.table(sigtab,"sigtab.txt",sep=",", row.names = TRUE)

#To manke the figures
DeSeq <- read.csv("DESeqResultsNew.csv", na.strings = c("","NA"), header=TRUE)
str(DeSeq)
DeSeq$High_low <- factor(DeSeq$High_low)
DeSeq$Species <- factor(DeSeq$Species)
levels(DeSeq$Species)
levels(DeSeq$Species) <- list("Bacteroides"="Bacteroides", "Bibersteinia"="Bibersteinia", "Helcococcus ovis"="Helcococcus ovis", "Moraxella (1)"="Moraxella (1)", "Moraxella (2)"="Moraxella (2)", "Moraxella boevrei DSM 14165 (1)"="Moraxella boevrei DSM 14165 (1)",
                              "Moraxellaceae (1)"="Moraxellaceae (1)","Mycoplasma (1)"="Mycoplasma (1)", "Mycoplasma (2)"="Mycoplasma (2)", "Mycoplasma alkalescens 14918"="Mycoplasma alkalescens 14918 ",
                              "Mycoplasma arginini"="Mycoplasma arginini", "Streptococcus"="Streptococcus" , "Trueperella pyogenes" ="Trueperella pyogenes",
                              "uncultured Parvimonas"="uncultured Parvimonas", "uncultured Bergeyella"="uncultured Bergeyella", "Clostridium sensu stricto 1"="Clostridium sensu stricto 1",
                              "Hydrogenophaga"="Hydrogenophaga", "Luteimonas"="Luteimonas", "Moraxella boevrei DSM 14165 (2)"="Moraxella boevrei DSM 14165 (2)", "Moraxellaceae (2)"="Moraxellaceae (2)",
                              "Mycoplasma bovirhinis"="Mycoplasma bovirhinis", "Salinicoccus"="Salinicoccus", "uncultured Gemmobacter"="uncultured Gemmobacter")

ggplot(data = DeSeq,aes(x = Species, y = log2FoldChange, group = factor(High_low))) + coord_flip() +
  geom_bar(stat = "identity", aes(fill = factor(High_low)), position = position_dodge(width = 0.9)) +
  labs(fill= "Diagnosis") +
  theme_bw()+
  ylab("Log2 Fold Change") +xlab ("Differentially ASVs") +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(8, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) 
  

#Functions
phyloseq_to_df <- function(physeq, addtax = T, addtot = F, addmaxrank = F, sorting = "abundance"){
  
  # require(phyloseq)
  
  ## Data validation
  if(any(addtax == TRUE || sorting == "taxonomy")){
    if(is.null(phyloseq::tax_table(physeq, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }
  }
  
  ## Prepare data frame
  if(taxa_are_rows(physeq) == TRUE){
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), phyloseq::otu_table(physeq), stringsAsFactors = F)
  } else {
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), t(phyloseq::otu_table(physeq)), stringsAsFactors = F)
  }
  
  ## Check if the sample names were silently corrected in the data.frame
  if(any(!phyloseq::sample_names(physeq) %in% colnames(res)[-1])){
    if(addtax == FALSE){
      warning("Warning: Sample names were converted to the syntactically valid column names in data.frame. See 'make.names'.\n")
    }
    
    if(addtax == TRUE){
      stop("Error: Sample names in 'physeq' could not be automatically converted to the syntactically valid column names in data.frame (see 'make.names'). Consider renaming with 'sample_names'.\n")
    }
  }
  
  ## Add taxonomy
  if(addtax == TRUE){
    
    ## Extract taxonomy table
    taxx <- as.data.frame(phyloseq::tax_table(physeq), stringsAsFactors = F)
    
    ## Reorder taxonomy table
    taxx <- taxx[match(x = res$OTU, table = rownames(taxx)), ]
    
    ## Add taxonomy table to the data
    res <- cbind(res, taxx)
    
    ## Add max tax rank column
    if(addmaxrank == TRUE){
      
      ## Determine the lowest level of taxonomic classification
      res$LowestTaxRank <- get_max_taxonomic_rank(taxx, return_rank_only = TRUE)
      
      ## Reorder columns (OTU name - Taxonomy - Max Rank - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), "LowestTaxRank", phyloseq::sample_names(physeq))]
      
    } else {
      ## Reorder columns (OTU name - Taxonomy - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), phyloseq::sample_names(physeq))]
      
    } # end of addmaxrank
  }   # end of addtax
  
  ## Reorder OTUs
  if(!is.null(sorting)){
    
    ## Sort by OTU abundance
    if(sorting == "abundance"){
      otus <- res[, which(colnames(res) %in% phyloseq::sample_names(physeq))]
      res <- res[order(rowSums(otus, na.rm = T), decreasing = T), ]
    }
    
    ## Sort by OTU taxonomy
    if(sorting == "taxonomy"){
      taxtbl <- as.data.frame( phyloseq::tax_table(physeq), stringsAsFactors = F )
      
      ## Reorder by all columns
      taxtbl <- taxtbl[do.call(order, taxtbl), ]
      # taxtbl <- data.table::setorderv(taxtbl, cols = colnames(taxtbl), na.last = T)
      res <- res[match(x = rownames(taxtbl), table = res$OTU), ]
    }
  }
  
  ## Add OTU total abundance
  if(addtot == TRUE){
    res$Total <- rowSums(res[, which(colnames(res) %in% phyloseq::sample_names(physeq))])
  }
  
  rownames(res) <- NULL
  return(res)
}


  

