library(afex)
library(lme4)
library(emmeans)
library(lubridate)
library(ggplot2)
library("cowplot")
theme_set(theme_grey())
#install.packages("jtools")
library(jtools)
library(ggpubr)
library(sjstats)

rm(list = ls ())

setwd("~/Desktop/eunice/MSc/Thesis/Scripts/PaperScript/Metadata/") 

metadata <- read.csv("16Smetadata.csv", na.strings = c("","NA"), header=TRUE)

#assign numerical values to factors
metadata$BRD <- as.factor(metadata$BRD)
levels(metadata$BRD) <- list("Healthy"="Healthy", "BRD"="BRD")
metadata$PenCode <- as.factor(metadata$PenCode)

str(metadata)

# Converting date of collection to numeric values
metadata$date <- metadata$Date.Collection
metadata$Date.Collection <- as.Date(metadata$Date.Collection, "%m/%d/%y")
d<- as.Date('12/31/2020', "%m/%d/%y") #use to calculate the days
metadata$Date.Collection <- as.Date(d) -as.Date(metadata$Date.Collection) 
metadata$Date.Collection <- as.numeric(metadata$Date.Collection)
str(metadata$Date.Collection)
#the highest day value is the date of the samples collected first 

plot(metadata$Date.Collection, metadata$Age)

## Dependent factors in the model
str(metadata)
#1. Observed OTUs
#2. Chao1 (measures richness of the environment)
#3. Pielou_e (measures evenness)
#4. Faith_pd (phylogenetic diversity)

set_sum_contrasts() # important for afex

# full model
str(metadata)
#For dependent variable Observed OTUs
#install.packages("piecewiseSEM")

M1 <- mixed(observed_otus ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1)

M2 <- mixed(pielou_e ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)

M3 <- mixed(chao1 ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)

M4 <- mixed(faith_pd ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4)

#checking assumptions
# way 1:
plot(M1$full_model)
plot(M2$full_model)
plot(M3$full_model)
plot(M4$full_model)

# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model))
qqnorm(residuals(M2$full_model))
qqnorm(residuals(M3$full_model))
qqnorm(residuals(M4$full_model))

# getting the means 
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_16S.1<- emmeans(M1, "Date.Collection")
emm_16S.1

# interpreting results
# BRD plots
a <- afex_plot(M1, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") +
  labs(y = "Observed ASVs", x = "Health Status") +
  theme(legend.position="none") 

c <- afex_plot(M2, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") +
  labs(y = "Evenness (Pielou)", x = "Health Status") +
  theme(legend.position="none")

b <- afex_plot(M3, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") +
  labs(y = "Chao 1", x = "Health Status") +
  theme(legend.position="none")

d <- afex_plot(M4, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") +
  labs(y = "Faith_pd", x = "Health Status") +
  theme(legend.position="none")


plot_grid(a, b, c , d, labels = "AUTO")


### Date collection plots
A <- ggplot(data = metadata, aes(x = Age, y = observed_otus)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("Observed ASVs") +xlab ("Animal Age (month)") +
  theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
  theme(text=element_text(family="Times New Roman"))

B <- ggplot(data = metadata, aes(x = Age, y = chao1)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("Chao1") +xlab ("Animal Age (month)") +
  theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
  theme(text=element_text(family="Times New Roman"))

C <- ggplot(data = metadata, aes(x = Age, y = pielou_e)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("Evenness (Pielou)") +xlab ("Animal Age (month)") +
  theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
  theme(text=element_text(family="Times New Roman"))

D <- ggplot(data = metadata, aes(x = Age, y = faith_pd)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("Faith_pd") +xlab ("Animal Age (month)") +
  theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
  theme(text=element_text(family="Times New Roman"))

ggarrange(A,B,C,D, labels = c("a", "b", "c", "d"),
          ncol = 2, nrow=2, font.label = list(size = 20))


e <- ggplot(data = metadata, aes(x = Date.Collection, y = observed_otus)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("Observed ASVs") +xlab ("Days relative to the end of the study") +
  theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
theme(text=element_text(family="Times New Roman"))

f <- ggplot(data = metadata, aes(x = Date.Collection, y = chao1)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("Chao1") +xlab ("Days relative to the end of the study") +
  theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
theme(text=element_text(family="Times New Roman"))

g <- ggplot(data = metadata, aes(x = Date.Collection, y = pielou_e)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("Evenness (Pielou)") +xlab ("Days relative to the end of the study") +
  theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
theme(text=element_text(family="Times New Roman"))

h <- ggplot(data = metadata, aes(x = Date.Collection, y = faith_pd)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("Faith_pd") +xlab ("Days relative to the end of the study") +
  theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
theme(text=element_text(family="Times New Roman"))

ggarrange(e,f,g,h, labels = c("a", "b", "c", "d"),
          ncol = 2, nrow=2, font.label = list(size = 20))

## Correlation between average temperature and date of collection

str(metadata)
cordata = metadata[,c(2,4)]
corr <- round(cor(cordata), 1)
corr

str(cordata)
cor(cordata$Date.Collection, cordata$Ave.temp)
cor.test(cordata$Date.Collection, cordata$Ave.temp)

library("ggpubr")
ggscatter(cordata, x = "Date.Collection", y = "Ave.temp", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Days relative to the end of the study", ylab = "Average temperature (°F)")


## testing temperature
str(meta)
set_sum_contrasts() 
M1 <- mixed(observed_otus ~ Ave.temp + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1)

M2 <- mixed(pielou_e ~ Ave.temp + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)

M3 <- mixed(chao1 ~ Ave.temp + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)

M4 <- mixed(faith_pd ~ Ave.temp + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4)

##Checking assumptions
#checking assumptions
# way 1:
plot(M1$full_model)
plot(M2$full_model)
plot(M3$full_model)
plot(M4$full_model)

# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model))
qqnorm(residuals(M2$full_model))
qqnorm(residuals(M3$full_model))
qqnorm(residuals(M4$full_model))

# Plots
i <- ggplot(data = metadata, aes(x = Ave.temp, y = observed_otus)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("Observed ASVs") +xlab ("Average daily temperature (°F)") +
  theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
theme(text=element_text(family="Times New Roman"))

j <- ggplot(data = metadata, aes(x = Ave.temp, y = chao1)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("Chao1") +xlab ("Average daily temperature (°F)") +
  theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
theme(text=element_text(family="Times New Roman"))

k <- ggplot(data = metadata, aes(x = Ave.temp, y = pielou_e)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("Evenness (Pielou)") +xlab ("Average daily temperature (°F)") +
  theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
theme(text=element_text(family="Times New Roman"))

l <- ggplot(data = metadata, aes(x = Ave.temp, y = faith_pd)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("Faith") +xlab ("Average daily temperature (°F)") +
  theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
theme(text=element_text(family="Times New Roman"))

ggarrange(i,j,k,l, labels = c("a", "b", "c", "d"),
          ncol = 2, nrow=2, font.label = list(size = 20))
