library(ggplot2)
library(vegan)
library(dplyr)
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
#install.packages("randomForest")
library(randomForest)
library(knitr)
library(qiime2R)
library(tidyr) #for separate function
library(naniar)# Ffor replace all function
library(ggpubr)

## Random Forest for 16S sequencing data
setwd("~/Desktop/eunice/MSc/Thesis/Qiime/RandomForest/124Samples/")
metadata <- read.csv("~/Desktop/eunice/MSc/Thesis/Scripts/PaperScript/Metadata_BRD/RandomF_124samples.csv", na.strings = c("","NA"), header=TRUE)
metadata2 <- metadata
str(metadata2)
metadata$BRD <- as.factor(metadata$BRD)
#metadata$Date.Collection <- mdy(metadata$Date.Collection)
order_groups <- metadata$ID
rownames(metadata) <- metadata[,1]
metadata = metadata[,-c(1)]
str(metadata)

ASVs <- read_qza("~/Desktop/eunice/MSc/Thesis/Qiime/Samples/OnlySamples/Filtered/Qiime/table-filtered2.qza")
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #18010 ASVs
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
ASV_table <- t(ASV_table)


##Adding taxonomy
#Taxonomy of each OTU
tax <- read_qza("~/Desktop/eunice/MSc/Thesis/Qiime/Samples/OnlySamples/Filtered/Qiime/taxonomy.qza")
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
tax.final$Family <- sub("D_9__*", "", tax.final[,4])
tax.final$Genus <- sub("D_0__*", "", tax.final[,5])
tax.final$Genus <- sub("D_1__*", "", tax.final[,5])
tax.final$Genus <- sub("D_2__*", "", tax.final[,5])
tax.final$Genus <- sub("D_3__*", "", tax.final[,5])
tax.final$Genus <- sub("D_4__*", "", tax.final[,5])
tax.final$Genus <- sub("D_5__*", "", tax.final[,5])
tax.final$Genus <- sub("D_9__*", "", tax.final[,5])
tax.final$Species <- sub("D_0__*", "", tax.final[,6])
tax.final$Species <- sub("D_1__*", "", tax.final[,6])
tax.final$Species <- sub("D_2__*", "", tax.final[,6])
tax.final$Species <- sub("D_3__*", "", tax.final[,6])
tax.final$Species <- sub("D_4__*", "", tax.final[,6])
tax.final$Species <- sub("D_5__*", "", tax.final[,6])
tax.final$Species <- sub("D_9__*", "", tax.final[,6])
#write.table(tax.final,"taxonomyNasal.txt",sep=",", row.names = FALSE) 
TaxASV <- merge(tax.final, ASVkey, by.x = 0, by.y = "ASVstring")
row.names(TaxASV) <- TaxASV[,10]
TaxASV = TaxASV[,-c(1,9:10)]
#write.table(TaxASV,"TaxASV.txt",sep=",", row.names = FALSE)

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
physeq_deseq ##17931 OTUS
## Random forest, we want to make the model so it classify samples based on the healthy status: BRD and Healthy
# Step 1: subset the samples based on healthy or sick

colnames(sample_data(physeq_deseq))

#Prune OTU that have an abundance <0.001
ntaxa(physeq_deseq) #total of 17931 ASVs

# Set prunescale 
prunescale = 0.0001
minlib = 40420 # rarefying to 1500 reads
# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean <- taxa_sums(physeq_deseq)/nsamples(physeq_deseq)
sites.prune <- prune_taxa(tax.mean > prunescale*minlib, physeq_deseq)

sites.prune
#1257 when working with 124 sampels
# Make a dataframe of training data with OTUs as column and samples as rows
predictors <- (otu_table(sites.prune))
dim(predictors) # we have 124 samples and 1257 OTUs 

# Make one column for our outcome/response variable 
response <- as.factor(sample_data(sites.prune)$BRD) ##classification
response

# Combine them into 1 data frame
rf.data <- data.frame(response, predictors) ## complete dataframe with the sample ID, classification and OTU table
str(rf.data)
#write.table(rf.data,"OtutableRF.txt",sep=",", row.names = TRUE) 

### other way
# Setting the testing and training set
str(metadata)
set.seed(123)
sample <- sample.int(n = nrow(metadata), size = floor(.60*nrow(metadata)), replace = F)
train <- metadata[sample, ]
str(train)
train %>% tally()
train %>% count(BRD)

test  <- metadata[-sample, ]
test %>% tally()
test %>% count(BRD)
str(test)
#write.table(test,"test.124.txt",sep=",", row.names = TRUE)

## performing random forest for the 16S data
train16S <- merge(train, rf.data, by.x = 0, by.y = 0)
rownames(train16S) <- train16S[,1]
train16S = train16S[,-c(1:3)]

S16 <- randomForest(response ~ .,data=train16S, importance=TRUE)
print(S16)
trainS16 <- as.data.frame(S16[["predicted"]])
trainS16$ID <- rownames(trainS16)
str(trainS16)
trainR <- merge(trainS16, metadata, by.x = "ID", by.y = 0)
#write.table(trainR,"trainR16S.124.txt",sep=",", row.names = TRUE) 

#testing data
test16S <- merge(test, rf.data, by.x = 0, by.y = 0)
rownames(test16S) <- test16S[,1]
test16S = test16S[,-c(1:3)]

predS16.2 = predict(S16, newdata=test16S)#predictions
predS16.2
pS16.2 <- as.data.frame(predS16.2)
pS16.2$ID <- rownames(pS16.2)
pS16.2 %>% count(predS16.2)
aS16 <- table(test16S$response, pS16.2$predS16.2)
abS16 <- as.data.frame(aS16)

#Confusion matrix
matrixS16 = aS16
matrixS16

#Misclasification rate
test16S %>%
  mutate(lda.pred = (pS16.2$pred)) %>%
  summarise(lda.error = mean(response != lda.pred))
rm(a)
a <- test16S %>%
  mutate(lda.pred = (pS16.2$pred)) %>%
  summarise(lda.error = mean(response != lda.pred))
View(a$lda.pred)
#Factor importance
varImpPlot(S16)
varImpPlot(S16, sort=TRUE,scale =FALSE, n.var=min(9,nrow(S16$importance)))

##Get the importance of the predictors in the randomForest model
imporatantASV <- as.data.frame(S16[["importance"]])
imporatantASV$ASV <- rownames(imporatantASV)
imporatantASV <- merge(imporatantASV, TaxASV, by.x = "ASV", by.y = 0)
str(imporatantASV)
imporASV <- subset(imporatantASV, MeanDecreaseAccuracy > 0.001)

str(test16S)
predS16 = predict(S16, newdata=test16S,type='prob')
predS16
pS16 <- as.data.frame(predS16)
pS16$ID <- rownames(pS16)
pS16
colnames(pS16) <- c("Prob.BRD16S", "Prob.Healthy16S", "ID")


###QPCR data
qPCR <- read.csv("~/Desktop/eunice/MSc/Thesis/Scripts/PaperScript/Metadata_BRD/qPCR_124samples.csv", na.strings = c("","NA"), header=TRUE)
str(qPCR)
qPCR$BRD <- as.factor(qPCR$BRD)
rownames(qPCR) <- qPCR[,1]
qPCR =qPCR[,-c(1)]
str(qPCR)

#RandomForest
trainqPCR <- merge(train, qPCR, by.x = 0, by.y = 0)
rownames(trainqPCR) <- trainqPCR[,1]
trainqPCR = trainqPCR[,-c(1:3)]

qpcr <- randomForest(BRD.y~ .,data=trainqPCR, importance=TRUE)
print(qpcr)
trainqpcr <- as.data.frame(qpcr[["predicted"]])
trainqpcr$ID <- rownames(trainqpcr)
str(trainqpcr)
trainqpcr <- merge(trainqpcr, metadata, by.x = 0, by.y = 0)
#write.table(trainR,"trainRqPCRnolog.124.txt",sep=",", row.names = TRUE) 

#testing set
testqPCR <- merge(test, qPCR, by.x = 0, by.y = 0)
rownames(testqPCR) <- testqPCR[,1]
testqPCR = testqPCR[,-c(1:3)]

predQPCR = predict(qpcr, newdata=testqPCR) #predictions
predQPCR
predQPCR <- as.data.frame(predQPCR)
predQPCR$ID <- rownames(predQPCR)
predQPCR %>% count(predQPCR)
aPCR <- table(testqPCR$BRD, predQPCR$predQPCR)
Apqpcr <- as.data.frame(aPCR)


#Confusion matrix
matrixqpcr = Apqpcr
matrixqpcr

#Misclasification rate
testqPCR %>%
  mutate(lda.pred = (predQPCR$pred)) %>%
  summarise(lda.error = mean(BRD.y != lda.pred))

#Factor importance
varImpPlot(qpcr)
varImpPlot(qpcr, sort=TRUE,scale =FALSE, n.var=min(10,nrow(qpcr$importance)))


##Get the importance of the predictors in the randomForest model
imporatantASV <- as.data.frame(qpcr[["importance"]])
#write.table(imporatantASV,"imporatantqPCR16S.txt",sep=",", row.names = TRUE)


predqpcr = predict(qpcr, newdata=testqPCR, type='prob') #predictions
predqpcr
pqpcr <- as.data.frame(predqpcr)
pqpcr$ID <- rownames(pqpcr)
pqpcr
colnames(pqpcr) <- c("Prob.BRD", "Prob.Healthy", "ID")


#plot decision Tree random forest
pqpcr
pred <- merge(pqpcr, metadata, by.x = "ID", by.y = 0)
pred <- merge(pred, pS16, by.x = "ID", by.y = "ID")
pred <- merge(pred, predQPCR, by.x = "ID", by.y = "ID")
pred <- merge(pred, pS16.2, by.x = "ID", by.y = "ID")
write.table(pred,"RandomForestpred2.txt",sep=",", row.names = TRUE) 

M3 <- read.csv("~/Desktop/eunice/MSc/Thesis/Qiime/RandomForest/124Samples/RandomForest.3M.csv", na.strings = c("","NA"), header=TRUE)
str(M3)
M3$classification <- as.factor(M3$classification)

#Prediction for healthy animals
str(pred)

ggplot() + 
  geom_point(data=pred, aes(x=Prob.BRD, y=Prob.BRD16S, size=3,color=BRD)) + 
  theme_classic() + 
  geom_point(data=M3, aes(x=Prob.BRD, y=Prob.BRD16S, size=3, shape=classification)) + 
  #geom_text(aes(label=ID),hjust=0.5, vjust=0) +
  guides(size=FALSE) +
  guides(fill=FALSE) +
  #geom_point(data= wuaxis, aes(x=groupX, y=groupX, shape=wuaxis4, size=5, fill=wuaxis4)) +
  scale_shape_manual(values=c(2, 8))+
  labs(color= "Visual Diagnosis") +
  labs(fill= "Health Status") +
  labs(shape= "Agreement 3 methods") +
  #geom_segment(data= wuaxis, aes(x=wuaxis1, y=wuaxis2, xend=groupX, yend=groupY, color= wuaxis4), size = .05) +
  labs(x='Predicted BRD-qPCR', y= 'Predicted BRD-ASV Table') +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) 

levels(test$BRD)
