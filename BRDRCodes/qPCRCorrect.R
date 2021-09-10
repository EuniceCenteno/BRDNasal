library(afex)
library(lme4)
library(emmeans)
library(lubridate)
library(psych)
library(sjstats)
library(tidyverse)
library(ggfortify)
library(pwr)
library(ggpubr)

# clear memory
#rm(list = ls ())

# load data
setwd("~/Desktop/eunice/MSc/Thesis/Scripts/PaperScript/Metadata/")

metadata <- read.csv("qPCRMetadata copy.csv", na.strings = c("","NA"), header=TRUE)
#assign numerical values to factors
metadata$BRD <- as.factor(metadata$BRD)
levels(metadata$BRD) <- list("Healthy"="Healthy", "BRD"="BRD")
metadata$PenCode <- as.factor(metadata$PenCode)
metadata$MbPre_Abs <- as.factor(metadata$MbPre_Abs)
metadata$MhPres_Abs <- as.factor(metadata$MhPres_Abs)
metadata$HsPres_Abs <- as.factor(metadata$HsPres_Abs)
metadata$PmPres_Abs <- as.factor(metadata$PmPres_Abs)

# Converting date to numeric value
metadata$date <- metadata$Date.Collection
metadata$Date.Collection <- as.Date(metadata$Date.Collection, "%m/%d/%y")
d<- as.Date('12/31/2020', "%m/%d/%y") #use to calculate the days
metadata$Date.Collection <- as.Date(d) -as.Date(metadata$Date.Collection) 
metadata$Date.Collection <- as.numeric(metadata$Date.Collection)
str(metadata$Date.Collection)
#the highest day value is the date of the samples collected first 

str(metadata)

#levels(metadata$Date)
plot(metadata$Date.Collection, metadata$Age)

##calcutating BRD-associated relative abundance based on the 16S counts
str(metadata)
metadata$Mb.Rel <- metadata$Mbovis_copies / metadata$X16S_copies
metadata$Pm.Rel <- metadata$Pm_copies / metadata$X16S_copies
metadata$Hs.Rel <- metadata$Hs_copies / metadata$X16S_copies
metadata$Mh.Rel <- metadata$Mh_copies / metadata$X16S_copies

## Dependent factors in the model
str(metadata)
## Dependent variables
#1. Mycoplasma bovis copies
#2. Pasteurella multocida copies
#3. Histophilus somni copies
#4. Mannheimia haemolytica copies
#5. 16S rRNA gene copies
#6. Mycoplasma bovis - relative abundance based on 16S rRNA gene copies
#7. Pasteurella multocida - relative abundance based on 16S rRNA gene copies
#8. Histophilus somni - relative abundance based on 16S rRNA gene copies
#9. Mannheimia haemolytica - relative abundance based on 16S rRNA gene copies
#10. Mycoplasma bovis - relative abundance based on the 4 bacteria total gene copies
#11. Pasteurella multocida - relative abundance based on the 4 bacteria total gene copies
#12. Histophilus somni - relative abundance based on the 4 bacteria total gene copies
#13. Mannheimia haemolytica - relative abundance based on the 4 bacteria total gene copies

set_sum_contrasts() # important for afex

# full model
str(metadata)

## 16S data
str(metadata)
S1 <- mixed(X16S_copies ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(S1)

#Checking assumptios
plot(S1$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(S1$full_model))

##Transform data
metadata <- mutate(metadata, X16S_log = log10(X16S_copies + 1))
metadata <- mutate(metadata, X16S_sqrt = sqrt(X16S_copies + 0.5))

S1.1 <- mixed(X16S_log ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(S1.1)
emm_S1.1<- emmeans(S1.1, "BRD")
emm_S1.1

S1.2 <- mixed(X16S_sqrt ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
              control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(S1.2)

#checking
plot(S1$full_model)
plot(S1.1$full_model)
plot(S1.2$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(S1$full_model))
qqnorm(residuals(S1.1$full_model)) #better
qqnorm(residuals(S1.2$full_model))

str(metadata)

afex_plot(S1.1, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 4), mapping="color") +
  theme_bw() + theme(legend.position="bottom") +
  labs(y = "16S rRNA gene copy (log10)", x = "Health Status") +
  theme(legend.position="none")  

#For dependent variable ---------qPCR
M1 <- mixed(Mbovis_copies ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1)

M2 <- mixed(Pm_copies ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)

M3 <- mixed(Hs_copies ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)

M4 <- mixed(Mh_copies ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4)

#checking assumptions
# way 1:
plot(M1$full_model)
plot(M2$full_model)
plot(M3$full_model)
plot(M4$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model)) #weird
qqnorm(residuals(M2$full_model)) #weird
qqnorm(residuals(M3$full_model))
qqnorm(residuals(M4$full_model)) #weird

##transformation
metadata <- mutate(metadata, Mbovis_log = log10(Mbovis_copies + 1))
metadata <- mutate(metadata, Mbovis_sqrt = sqrt(Mbovis_copies + 0.5))

metadata <- mutate(metadata, Pm_log = log10(Pm_copies + 1))
metadata <- mutate(metadata, Pm_sqrt = sqrt(Pm_copies + 0.5))

metadata <- mutate(metadata, Mh_log = log10(Mh_copies + 1))
metadata <- mutate(metadata, Mh_sqrt = sqrt(Mh_copies + 0.5))

## Test new variables
#For dependent variable 
M1.1 <- mixed(Mbovis_log ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1.1)

M1.2 <- mixed(Mbovis_sqrt ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1.2)

M2.1 <- mixed(Pm_log ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2.1)

M2.2 <- mixed(Pm_sqrt ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
              control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2.2)

M4.1 <- mixed(Mh_log ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
              control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4.1)

M4.2 <- mixed(Mh_sqrt ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
              control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4.2)

#checking assumptions
# way 1:
plot(M1$full_model)
plot(M1.1$full_model) #better
plot(M1.2$full_model)

plot(M2$full_model)
plot(M2.1$full_model)#better
plot(M2.2$full_model)  

plot(M4$full_model)
plot(M4.1$full_model) #better
plot(M4.2$full_model)  

# this is for testing the normality of the residuals
qqnorm(residuals(M1$full_model)) 
qqnorm(residuals(M1.1$full_model)) #better log
qqnorm(residuals(M1.2$full_model))

qqnorm(residuals(M2$full_model)) 
qqnorm(residuals(M2.1$full_model)) #better log
qqnorm(residuals(M2.2$full_model)) 

qqnorm(residuals(M4$full_model)) 
qqnorm(residuals(M4.1$full_model)) #better log
qqnorm(residuals(M4.2$full_model)) 

## new variables, Mbovis_log, Pm_sqrt, and Mh_log
M1.1 <- mixed(Mbovis_log ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
              control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1.1)
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_M1.1<- emmeans(M1.1, "BRD")
emm_M1.1

M2.1 <- mixed(Pm_log ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
              control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2.1)

M3 <- mixed(Hs_copies ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_M3<- emmeans(M3, "BRD")
emm_M3

M4.1 <- mixed(Mh_log ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
              control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4.1)
emm_M4.1<- emmeans(M4.1, "BRD")
emm_M4.1


# interpreting results
afex_plot(M1.1, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 6), mapping="color") +
  theme_bw(base_size = 15) + 
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +  
  labs(y = "M. bovis copies (log10)", x = "Health Status") +
  theme(legend.position="none") 

afex_plot(M3, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 6), mapping="color") +
  theme_bw(base_size = 15) + 
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) + 
  labs(y = "H. somni copies", x = "Health Status") +
  theme(legend.position="none")

afex_plot(M4.1, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 6), mapping="color") +
  theme_bw(base_size = 15) + 
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(y = "M. haemolytica copies (log10)", x = "Health Status") +
  theme(legend.position="none")

## Relative abundance values
str(metadata)
m1 <- mixed(Mb.Rel ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(m1)

m2 <- mixed(Pm.Rel ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(m2)
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_m2<- emmeans(m2, "BRD")
emm_m2

m3 <- mixed(Hs.Rel ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(m3)

m4 <- mixed(Mh.Rel ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(m4)

###Plots
afex_plot(m2, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 6), mapping="color") +
  theme_bw(base_size = 15) + 
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) + 
  labs(y = "P. multocida- Rel.Abundance", x = "Health Status") +
  theme(legend.position="none")


