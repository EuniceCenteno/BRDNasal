BETA diversity Swabs

library(vegan) 
library(ggplot2)
library(ggpubr)
library(data.table)
library(phyloseq)
library(qiime2R)
library(tidyr)
library(naniar)
library(raster)
#transpose function

setwd("~/Desktop/eunice/MSc/Thesis/Qiime/Samples/OnlySamples/Filtered/Qiime/") #sets new working directory for Windows systems (remember to replace … with your filepath)
metadata <-read.csv("16SBetametadata.csv", na.strings = c("","NA"), header=TRUE)
#OTU table (shared file)
#The OTU table as exported from qiime has a pound sign before the header row. You need to delete that pound sign in a text editor.
str(metadata)
metadata$BRD <- factor(metadata$BRD) 
metadata$PenCode <- factor(metadata$PenCode)
order_groups <- metadata$ID
row.names(metadata) = metadata[,1]
metadata = metadata[,-1]

#OTU table, we use the rarified table
ASVs <- read_qza("rarefied_table.qza")
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
tax.clean$OTUs <- rownames(tax.clean)
#Would be good to check here to make sure the order of the two data frames was the same. You should do this on your own.

###Remove all the OTUs that don't occur in our OTU.clean data set
tax.final = tax.clean[row.names(tax.clean) %in% row.names(ASV_s),]

tax.final$Phylum <- sub("D_1__*", "", tax.final[,1])
tax.final$Class <- sub("D_2__*", "", tax.final[,2])
tax.final$Order <- sub("D_3__*", "", tax.final[,3])
tax.final$Family <- sub("D_9__*", "", tax.final[,4])
tax.final$Genus <- sub("D_10__*", "", tax.final[,5])
tax.final$Species <- sub("D_10__*", "", tax.final[,6])
#write.table(tax.final,"taxonomyNasal.txt",sep=",", row.names = FALSE) 
TaxASV <- merge(tax.final, ASVkey, by.x = 0, by.y = "ASVstring")
row.names(TaxASV) <- TaxASV[,10]
TaxASV = TaxASV[,-c(1,10)]
#write.table(TaxASV,"TaxASV.txt",sep=",", row.names = FALSE)

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
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]
#write.table(NewTax,"NewTax.txt",sep=",", row.names = TRUE)
NewASVTable <- prunetable
NewASVTable <- NewASVTable[,-c(2:9)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]

## Calculating distances based on Bray-curtis and Weighted Unifrac using the physeq object
dist.bray <- phyloseq::distance(physeq_deseq, method = "bray")
dist.bray <- as.dist(dist.bray)

## PERMANOVA
#Bray curtis
BC <- adonis(dist.bray ~ metadata$BRD, permutations = 999)
BC

## Weighted Unifrac
# We need to create first a tree using OTU and taxa table-- we do this by creating a phyloseq object 

OTU.physeq = otu_table(as.matrix(NewASVTable), taxa_are_rows=TRUE)
tax.physeq = tax_table(as.matrix(NewTax))

#We then merge these into an object of class phyloseq.
physeq = phyloseq(OTU.physeq, tax.physeq)
physeq

##Making the tree
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
##Merging the tree with the phyloseq object
physeq1 = merge_phyloseq(physeq, meta.physeq, random_tree)
physeq1 

OTU_Unifrac <- UniFrac(physeq1, weighted = TRUE) ## Calculating the weighted unifrac distances

## Permanova 
M3 <- adonis(OTU_Unifrac~ metadata$BRD, permutations = 999)
M3

## Dispersion test 
BRD_Wu <- betadisper(OTU_Unifrac, type = c("centroid"), group = metadata$BRD)
BRD_Wu
boxplot(BRD_Wu)
TukeyHSD(BRD_Wu) ## calculate the distance from the centroids to one group to another
plot(BRD_Wu)

#install.packages("usedist")
library(usedist)

pBRD_WU<- permutest(BRD_Wu, permutations = 999)
pBRD_WU# not significant 

#Bray_curtis dispersion
BRD_BC <- betadisper(dist.bray, type = c("centroid"), group = metadata$BRD)
BRD_BC
boxplot(BRD_BC)

pBRD_BC<- permutest(BRD_BC, permutations = 999)
pBRD_BC# not significant 

##O Make the plots with ellipses
## Calculating the centroids of the data based on Bray-curtis and weighted unifrac
#D33 centroids for treatment and season

## Weighted unifrac 
ordu.wt.uni <- ordinate(physeq1 , "PCoA", "unifrac", weighted=T)
wt.unifrac <- plot_ordination(physeq1, 
                              ordu.wt.uni, color="BRD") 
wuaxis1 <- wt.unifrac[["data"]][["Axis.1"]]
wuaxis1 <- as.data.frame(wuaxis1)
wuaxis1$number <- rownames(wuaxis1)
wuaxis2 <- wt.unifrac[["data"]][["Axis.2"]]
wuaxis2 <- as.data.frame(wuaxis2)
wuaxis2$number <- rownames(wuaxis2)
wuaxis3 <- wt.unifrac[["data"]][["Cattle.ID"]]
wuaxis3 <- as.data.frame(wuaxis3)
wuaxis3$number <- rownames(wuaxis3)
wuaxis4 <- wt.unifrac[["data"]][["BRD"]]
wuaxis4 <- as.data.frame(wuaxis4)
wuaxis4$number <- rownames(wuaxis4)

wuaxis <- merge(wuaxis3, wuaxis4, by.x = "number", by.y = "number")
wuaxis <- merge(wuaxis, wuaxis1, by.x = "number", by.y = "number")
wuaxis <- merge(wuaxis, wuaxis2, by.x = "number", by.y = "number")


wt.unifrac <- wt.unifrac + ggtitle("Weighted UniFrac") + geom_point(size = 2)
wt.unifrac <- wt.unifrac + theme_classic() 
print(wt.unifrac + stat_ellipse())

##other plot
str(wuaxis)

centroids <- aggregate(wuaxis[,4:5], list(Group=wuaxis$wuaxis4), mean)
colnames(centroids) <- c('BRD','groupX', 'groupY')

wuaxis <- merge(wuaxis, centroids, by.x = "wuaxis4", by.y = "BRD")

str(wuaxis)
a <- ggplot(wuaxis, aes(x=wuaxis1, y=wuaxis2, color=wuaxis4)) + 
  geom_point(size=2) + 
  theme_classic() + stat_ellipse() +
  guides(size=FALSE) +
  geom_point(data= wuaxis, aes(x=groupX, y=groupX, shape=wuaxis4, size=5, fill=wuaxis4)) +
  scale_shape_manual(values=c(25, 22))+
  labs(color= "Health Status") +
  labs(fill= "Health Status") +
  labs(shape= "Centroids") +
  #geom_segment(data= wuaxis, aes(x=wuaxis1, y=wuaxis2, xend=groupX, yend=groupY, color= wuaxis4), size = .05) +
  labs(x='Axis 1 (21.7%)', y= 'Axis 2  (13.4%)', caption = paste('Distance between centroids: 0.057')) +
  ggtitle("Weighted UniFrac") +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) 


#distance betweeen centroids
pointDistance(c(0.03168905, -0.0010942445), c(-0.02609687, 0.0009011425), lonlat=FALSE)

## Bray_curtis
ordu.bc <- ordinate(physeq1, "PCoA", "bray")
Bray <- plot_ordination(physeq1, 
                              ordu.bc, color="BRD") 
Bray1 <- Bray[["data"]][["Axis.1"]]
Bray1 <- as.data.frame(Bray1)
Bray1$number <- rownames(Bray1)
Bray2 <- Bray[["data"]][["Axis.2"]]
Bray2 <- as.data.frame(Bray2)
Bray2$number <- rownames(Bray2)
Bray3 <- Bray[["data"]][["Cattle.ID"]]
Bray3 <- as.data.frame(Bray3)
Bray3$number <- rownames(Bray3)
Bray4 <- Bray[["data"]][["BRD"]]
Bray4 <- as.data.frame(Bray4)
Bray4$number <- rownames(Bray4)

bray <- merge(Bray3, Bray4, by.x = "number", by.y = "number")
bray <- merge(bray, Bray1, by.x = "number", by.y = "number")
bray <- merge(bray, Bray2, by.x = "number", by.y = "number")

Bray <- Bray + ggtitle("Bray Curtis") + geom_point(size = 2)
Bray <- Bray + theme_classic()
print(Bray + stat_ellipse())

str(bray)

centroids2 <- aggregate(bray[,4:5], list(Group=bray$Bray4), mean)
colnames(centroids2) <- c('BRD','groupX', 'groupY')

bray <- merge(bray, centroids2, by.x = "Bray4", by.y = "BRD")
## Distance between centroids
pointDistance(c(0.03947763, -0.02950140), c(-0.03251099, 0.02429527), lonlat=FALSE)


b <- ggplot(bray, aes(x=Bray1, y=Bray2, color=Bray4)) + 
  geom_point(size=2) + 
  theme_classic() + stat_ellipse() +
  guides(size=FALSE) +
  geom_point(data= bray, aes(x=groupX, y=groupX, shape=Bray4, size=5, fill=Bray4)) +
  scale_shape_manual(values=c(25, 22))+
  labs(color= "Health Status") +
  labs(fill= "Health Status") +
  labs(shape= "Centroids") +
  labs(x='Axis 1 (16.2%)', y= 'Axis 2 (13.6%)', caption = paste('Distance between centroids: 0.089')) +
  ggtitle("Bray-Curtis dissimilarity") +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) 


ggarrange(a,b,labels = c("a", "b"),
          ncol = 2, font.label = list(size = 20))

#Function
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





