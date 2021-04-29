# This is a demo for running the co-occurrence analysis much, much faster

#make sure you have these libraries
#install.packages("zoo")
library(Hmisc)
library(plyr)
library(fdrtool)
library(intergraph)
library(ggplot2)
library(tidyr)
library(lubridate)
library(cooccur)
library(qiime2R)
library(phyloseq)
data(finches)
library(naniar) ##for replace_with_na_all function


#install.packages("remotes")
#remotes::install_github("jbisanz/qiime2R")

# this is the data
setwd("~/Desktop/eunice/Thesis/Qiime/Samples/OnlySamples/Filtered/network/")

#OTU table (shared file)
#The OTU table as exported from qiime has a pound sign before the header row. You need to delete that pound sign in a text editor.
metadata <- read.delim("~/Desktop/eunice/Thesis/Qiime/Samples/OnlySamples/Filtered/DESeq/DESeqmetadata.txt", sep = "\t", header = T, quote = "", stringsAsFactors = F)
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
# Step 1: subset the samples based on healthy or sick

#Pruning the data
# Set prunescale 
prunescale = 0.0001
minlib = 40420 # rarefying to 1500 reads
# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean <- taxa_sums(physeq_deseq)/nsamples(physeq_deseq)
sites.prune <- prune_taxa(tax.mean > prunescale*minlib, physeq_deseq)

sites.prune

prunetable4<- phyloseq_to_df(sites.prune, addtax = T, addtot = F, addmaxrank = F,
               sorting = "abundance")

#Subseting the OTU table based on healthy and sick
OTUtable <- prunetable4[,-c(2:9)]
row.names(OTUtable) <- OTUtable[,1]
OTUtable<- OTUtable[,-c(1)]
OTUtable <- t(OTUtable)

# Taxonomy table
TaxaPrune<- prunetable4[,-c(10:140)]
## Merging metadata en OTU table to then divide the data
str(metadata)
data <- merge(metadata, OTUtable, by.x = 0, by.y = 0)
str(data)

## Subset the data into healthy and sick
dataHealthy <- subset(data, BRD=="0")
dataHealthy<- dataHealthy[,-c(2:3)] #to remove any unnecessary columns
row.names(dataHealthy) <-dataHealthy[,1]
dataHealthy<- dataHealthy[,-c(1)]
dataHealthy <- t(dataHealthy)
#write.table(dataHealthy,"dataHealthy.txt",sep=",", row.names = TRUE) 
dataHealthy[dataHealthy > 1] <- 1 ## we transpose because for the co-occur function, it needs the rows be ASVs and colums sites (ID)
#write.table(dataHealthy,"dataHealthyPresence.txt",sep=",", row.names = TRUE) 

dataBRD <- subset(data, BRD=="1")
dataBRD<- dataBRD[,-c(2:3)] #to remove any unnecessary columns
row.names(dataBRD) <-dataBRD[,1]
dataBRD<- dataBRD[,-c(1)]
dataBRD<- t(dataBRD)
#write.table(dataBRD,"dataBRD.txt",sep=",", row.names = TRUE) 
dataBRD[dataBRD > 1] <- 1## we transpose because for the co-occur function, it needs the rows be ASVs and colums sites (ID)
#write.table(dataBRD,"dataBRDPresence.txt",sep=",", row.names = TRUE) 

## Running co-occur 
B <- cooccur(mat=dataBRD,
                           type="spp_site", thresh=TRUE,
                           spp_names=TRUE)
summary(H1co)
plot(H1co)

H <- cooccur(mat=dataHealthy,
             type="spp_site", thresh=TRUE,
             spp_names=TRUE)
summary(H1co)
plot(H1co) 


##Cooocurrence data Healthy
CooccurH <- read.csv("printHealthy.txt", na.strings = c("","NA"), header=TRUE)
CoocpronH <- read.csv("prob.tableHealthy.txt", na.strings = c("","NA"), header=TRUE)
### CooccurH contains only the positive and negative pairs
# p_lt =1 positive 
# p_gt= negative

positivepairs <-subset(CooccurH, p_lt > 0.9) 
positivepairs <- merge(TaxaPrune, positivepairs, by.x = "OTU", by.y = "sp1_name")
colnames(positivepairs) <- c("sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                      "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "sp2_name")
positivepairs <- merge(TaxaPrune,positivepairs, by.x = "OTU", by.y = "sp2_name")
colnames(positivepairs) <- c("sp2_name", "Phylum_2", "Class_2", "Order_2", "Family_2", "Genus_2", "Species_2", "Confidence_2", "OTUs_2",
                      "sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                      "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt")
#write.table(positivepairs,"positivepairsHealthy.txt",sep=",", row.names = FALSE) 
positiveH<-subset(CooccurH, p_lt > 0.999999) ##only the ones with high positive probability
positiveH <- subset(CooccurH, obs_cooccur > 45) ### probability that the species will be in at least 60% of the total samples
positiveH <- subset(CooccurH, prob_cooccur > 0.9) ### probability that the two species will be in the same site

str(positiveH)
PosTax <- merge(TaxaPrune, positiveH, by.x = "OTU", by.y = "sp1_name")
colnames(PosTax) <- c("sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                      "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "sp2_name")
PosTax <- merge(TaxaPrune,PosTax, by.x = "OTU", by.y = "sp2_name")
colnames(PosTax) <- c("sp2_name", "Phylum_2", "Class_2", "Order_2", "Family_2", "Genus_2", "Species_2", "Confidence_2", "OTUs_2",
                      "sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                      "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt")
#write.table(PosTax,"PosTax.txt",sep=",", row.names = FALSE) 

negativepairs <-subset(CooccurH, p_lt < 0.01)
negativepairs <-subset(negativepairs, obs_cooccur <1)
negativepairs <-subset(negativepairs, prob_cooccur <0.05)
negativepairs <- merge(TaxaPrune, negativepairs, by.x = "OTU", by.y = "sp1_name")
colnames(negativepairs) <- c("sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                             "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "sp2_name")
negativepairs <- merge(TaxaPrune,negativepairs, by.x = "OTU", by.y = "sp2_name")
colnames(negativepairs) <- c("sp2_name", "Phylum_2", "Class_2", "Order_2", "Family_2", "Genus_2", "Species_2", "Confidence_2", "OTUs_2",
                             "sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                             "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt")
write.table(negativepairs,"negativepairsHealthy.txt",sep=",", row.names = FALSE) 

negaH<-subset(CooccurH, p_lt < 0.05) ##only the ones with high positive probability
negaH <- merge(TaxaPrune, negaH, by.x = "OTU", by.y = "sp1_name")
colnames(negaH) <- c("sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                             "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "sp2_name")
negaH <- merge(TaxaPrune,negaH, by.x = "OTU", by.y = "sp2_name")
colnames(negaH) <- c("sp2_name", "Phylum_2", "Class_2", "Order_2", "Family_2", "Genus_2", "Species_2", "Confidence_2", "OTUs_2",
                             "sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                             "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt")
write.table(negaH,"negativeHealthy.txt",sep=",", row.names = FALSE) 

##Cooocurrence data BRD
CooccurB <- read.csv("printBRD.txt", na.strings = c("","NA"), header=TRUE)
CoocpronB <- read.csv("prob.tableBRD.txt", na.strings = c("","NA"), header=TRUE)
### CooccurH contains only the positive and negative pairs
# p_lt =1 positive 
# p_gt= negative

positivepairsB <-subset(CooccurB, p_lt > 0.9) 
positivepairsB <- merge(TaxaPrune, positivepairsB, by.x = "OTU", by.y = "sp1_name")
colnames(positivepairsB) <- c("sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                             "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "sp2_name")
positivepairsB <- merge(TaxaPrune,positivepairsB, by.x = "OTU", by.y = "sp2_name")
colnames(positivepairsB) <- c("sp2_name", "Phylum_2", "Class_2", "Order_2", "Family_2", "Genus_2", "Species_2", "Confidence_2", "OTUs_2",
                             "sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                             "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt")

#write.table(positivepairsB,"positivepairsB.txt",sep=",", row.names = TRUE) 
positiveB<-subset(CooccurB, p_lt > 0.999999) ##only the ones with high positive probability
positiveB <- subset(CooccurB, obs_cooccur > 34) ### probability that the species will be in at least 60% of the total samples
positiveB <- subset(CooccurB, prob_cooccur > 0.9) ### probability that the two species will be in the same site

str(positiveB)
PosTaxB <- merge(TaxaPrune, positiveB, by.x = "OTU", by.y = "sp1_name")
colnames(PosTaxB) <- c("sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                      "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "sp2_name")
PosTaxB <- merge(TaxaPrune,PosTaxB, by.x = "OTU", by.y = "sp2_name")
colnames(PosTaxB) <- c("sp2_name", "Phylum_2", "Class_2", "Order_2", "Family_2", "Genus_2", "Species_2", "Confidence_2", "OTUs_2",
                      "sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                      "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt")
#write.table(PosTaxB,"PosTaxB.txt",sep=",", row.names = TRUE) 

negativepairsB <-subset(CooccurB, p_lt < 0.01)
negativepairsB <-subset(negativepairsB, obs_cooccur <1)
negativepairsB <-subset(negativepairsB, prob_cooccur <0.05)
negativepairsB <- merge(TaxaPrune, negativepairsB, by.x = "OTU", by.y = "sp1_name")
colnames(negativepairsB) <- c("sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                             "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "sp2_name")
negativepairsB <- merge(TaxaPrune,negativepairsB, by.x = "OTU", by.y = "sp2_name")
colnames(negativepairsB) <- c("sp2_name", "Phylum_2", "Class_2", "Order_2", "Family_2", "Genus_2", "Species_2", "Confidence_2", "OTUs_2",
                             "sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                             "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt")
write.table(negativepairsB,"negativepairsBRD.txt",sep=",", row.names = FALSE) 

negaB<-subset(CooccurB, p_lt < 0.05) ##only the ones with high positive probability
negaB <- merge(TaxaPrune, negaB, by.x = "OTU", by.y = "sp1_name")
colnames(negaB) <- c("sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                     "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt", "sp2_name")
negaB <- merge(TaxaPrune,negaB, by.x = "OTU", by.y = "sp2_name")
colnames(negaB) <- c("sp2_name", "Phylum_2", "Class_2", "Order_2", "Family_2", "Genus_2", "Species_2", "Confidence_2", "OTUs_2",
                     "sp1_name", "Phylum_1", "Class_1", "Order_1", "Family_1", "Genus_1", "Species_1", "Confidence_1", "OTUs_1", "sp1", "sp2", "sp1_inc", "sp2_inc", 
                     "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", "p_gt")
write.table(negaB,"negativeBRD.txt",sep=",", row.names = FALSE) 


##Function
https://rdrr.io/github/vmikk/metagMisc/src/R/phyloseq_to_df.R

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
