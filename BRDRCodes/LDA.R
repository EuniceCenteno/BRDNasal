Discriminant Function Analysis

library(tidyverse)
library(modelr)
library(broom)
library(ISLR)
library(ROCR)
library(MASS)

#install.packages("ISLR")
## Test the data
setwd("~/Desktop/eunice/Thesis/qPCR/LDA_with0/")
qPCR <- read.csv("LDAMetadata.csv", na.strings = c("","NA"), header=TRUE) #data with 0
qPCR <- read.csv("LDAMetadata_no0.csv", na.strings = c("","NA"), header=TRUE) # data of the animals that tested positive for the 4 bacteria, no 0 included
str(qPCR)

##transfor the data to log10
qPCR <- mutate(qPCR, Mbovis_copies = log10(Mbovis_copies + 1))
qPCR <- mutate(qPCR, Pm_copies = log10(Pm_copies + 1))
qPCR <- mutate(qPCR, Hs_copies = log10(Hs_copies + 1))
qPCR <- mutate(qPCR, Mh_copies = log10(Mh_copies + 1))
qPCR <- mutate(qPCR, X16S_copies = log10(X16S_copies + 1))
QPCR<- qPCR

str(QPCR)
#change the columns depending on the variables you will use to construct the model
qPCR <- QPCR[c(1,7,8)]
str(qPCR)
set.seed(123)
#training_sample2 <- sample(c(TRUE, FALSE), nrow(qPCR), replace = T, prob = c(0.6,0.4))
sample <- sample.int(n = nrow(qPCR), size = floor(.60*nrow(qPCR)), replace = F)
train <- qPCR[sample, ]
train$rownames <- rownames(train)

test  <- qPCR[-sample, ]


#Apply LDA to the training set
str(train)
train2<- train[c(1:3)]
lda.BRD <- lda(BRD ~ ., train2)
lda.BRD #show results

#Only on LD1 axis
plot(lda.BRD)


# See if the model fits the data -- we use the predict function
#we have to use the vector we created when we run the LDA 
lda.BRD2 <- predict(lda.BRD, newdata=test)


#Validation, how well the model predicts the true positives and true negatives (%)
str(test)
#chao and shannon are number
test %>% tally()
test %>% count(BRD)
train2 %>% count(BRD)

Tlda <- lda.BRD2[["class"]]
Tlda <- as.data.frame(Tlda)
Tlda %>% tally()
Tlda %>% count(Tlda)

lda.cm <- table(test$BRD, lda.BRD2$class)
cm <- as.data.frame(lda.cm)
cm
#Clasification rate
LDA_model = lda.cm %>% prop.table() %>% round(3)
LDA_model ### percent of true positives and true negatives


#Misclasification rate
test %>%
  mutate(lda.pred = (lda.BRD2$class)) %>%
  summarise(lda.error = mean(BRD != lda.pred))
        

# Confusion matrices
#Confusion matrix
LDA_model = lda.cm
LDA_model

