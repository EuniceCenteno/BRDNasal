library(ggplot2)
library(tidyr) #separate function
library(reshape2) #melt function
library(dplyr)
library(naniar) # for replace_with_na_all function
library(data.table)
library(qiime2R)

setwd("~/Desktop/eunice/MSc/Thesis/Qiime/Controls/Mock/")

ASVs <- read_qza("table.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #54 ASVs

#####################################################################
######################################################################

##Adding taxonomy
#Taxonomy of each OTU
tax = read.table("taxonomy.tsv", header=TRUE, sep="\t")
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
#write.table(tax.final,"taxonomyMock2.txt",sep=",", row.names = FALSE) 


##Calculating abundance of the ASVs
ASV_table = t(ASV_table)
otu.summary <- prop.table(as.matrix(ASV_table), 1) 
str(otu.summary)
otu_abund <- colSums(otu.summary)
otu_abund2 <- as.data.frame(otu_abund)
otu.summary <- rbind(otu_abund, otu.summary)
str(otu.summary)
otu.summary_sorted <- otu.summary[,order(otu.summary[1,], decreasing = TRUE)]
str(otu.summary_sorted)
melt_otu <- reshape2::melt(otu.summary_sorted[, c(1:54)]) ###TOTAL NUMBER OF OTUS
str(melt_otu)
colnames(melt_otu) <- c("Sample", "ASV", "Abundance")
str(melt_otu)
levels(melt_otu$Sample)
#write.table(melt_otu,"ASVAbundanceMock.txt",sep=",", row.names = FALSE)

## Merging contaminates AVSs with the AVS from the empty swabs -- see if they are contaminated with it as well

##- Empty swabs Taxonomy
setwd("~/Desktop/eunice/MSc/Thesis/Qiime/Controls/Empty/")

ASVs <- read_qza("table.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #54 ASVs

#####################################################################
######################################################################

##Adding taxonomy
#Taxonomy of each OTU
tax = read.table("taxonomy.tsv", header=TRUE, sep="\t")
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
#write.table(tax.final,"taxonomyEmpty.txt",sep=",", row.names = FALSE)

## Combine the two taxa tables so see which ones are shared and unique
setwd("~/Desktop/eunice/MSc/Thesis/Qiime/Controls/Mock/")
Cont_AVS = read.table("ContaminantsMock.txt", header=TRUE, sep="\t")

setwd("~/Desktop/eunice/MSc/Thesis/Qiime/Controls/Empty/")
TaxaEmpty <- read.table("taxonomyEmpty.txt", header=TRUE, sep="\t")
TaxCont_ASV <- merge(Cont_AVS, TaxaEmpty, by.x = "OTUs", by.y = "OTUs")
## No shared ASVs, no contamination between Mock and Empty Swabs

### Now getting the samples taxonomy
setwd("~/Desktop/eunice/MSc/Thesis/Qiime/Samples/OnlySamples/Filtered/Qiime/")

ASVs <- read_qza("table-filtered2.qza.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #54 ASVs

#####################################################################
######################################################################

##Adding taxonomy
#Taxonomy of each OTU


tax <- read_qza("taxonomy.qza")
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
tax.clean$OTUs <- rownames(tax.clean)C
#Would be good to check here to make sure the order of the two data frames was the same. You should do this on your own.

###Remove all the OTUs that don't occur in our OTU.clean data set
tax.final = tax.clean[row.names(tax.clean) %in% row.names(ASV_s),]
#write.table(tax.final,"taxonomySamples.txt",sep=",", row.names = FALSE)

#merging the swabs and the samples taxonomy
EmptySwabs_ASV <- merge(TaxaEmpty, tax.final, by.x = "OTUs", by.y = "OTUs")
#write.table(EmptySwabs_ASV,"SharedASVEmpty_Samples.txt",sep=",", row.names = FALSE)

# Now comparing samples and Mock contaminants
MockSwabs_ASV <- merge(Cont_AVS, tax.final, by.x = "OTUs", by.y = "OTUs")
#write.table(MockSwabs_ASV,"SharedASVMock_Samples.txt",sep=",", row.names = FALSE)
 
#Counting how many shared ASVs are between then ASV and empty tubes and samples
setwd("~/Desktop/eunice/MSc/Thesis/Qiime/Samples/OnlySamples/NoFilter/")
shared = read.csv("SharedASVEmpty_Samples.csv", header=TRUE)
str(shared)
shared$Family <- factor(shared$Family)
shared$Phylum <- factor(shared$Phylum)
shared$Genus <- factor(shared$Genus)

shared %>% tally()
str(shared)
shared %>% count(Family, sort = TRUE)
a <- shared %>% count(Genus, sort = TRUE)

perG <- shared %>%
  group_by(Genus) %>%
  summarise(count = n() ) %>%
  mutate( prop = (count / sum(count))*100 )

sum(perG$count)
#299 ASVs
sum(perG$prop)
#100 

perF <- shared %>%
  group_by(Family) %>%
  summarise(count = n() ) %>%
  mutate( prop = (count / sum(count))*100 )

sum(perF$count)
#299 ASVs
sum(perF$prop)

