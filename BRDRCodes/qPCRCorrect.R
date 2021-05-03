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
setwd("~/Desktop/eunice/Thesis/Scripts/PaperScript/Metadata/")

metadata <- read.csv("qPCRMetadata.csv", na.strings = c("","NA"), header=TRUE)
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
plot(M2.1$full_model)
plot(M2.2$full_model) #better 

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
  labs(y = "M. bovis copies (log10)- All Samples", x = "Health Status") +
  theme(legend.position="none") 

afex_plot(M3, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 6), mapping="color") +
  theme_bw(base_size = 15) + 
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) + 
  labs(y = "H. somni copies- All Samples", x = "Health Status") +
  theme(legend.position="none")

afex_plot(M4.1, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 6), mapping="color") +
  theme_bw(base_size = 15) + 
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(y = "M. haemolytica copies (log10)-All Samples", x = "Health Status") +
  theme(legend.position="none")

## Relative abundance values
str(metadata)
m1 <- mixed(Mbovis_RelAbunBacteria ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(m1)
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'

m2 <- mixed(Pm_RelAbunBacteria ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(m2)


m3 <- mixed(Hs_RelAbunBacteria ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(m3)

m4 <- mixed(Mh_RelAbunBacteria ~ BRD + Date.Collection + Age + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(m4)

###Plots
afex_plot(m1, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 6), mapping="color") +
  theme_bw(base_size = 15) + 
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +  
  labs(y = "M. bovis- Rel.Abun (All Samples)", x = "Health Status") +
  theme(legend.position="none") 

afex_plot(m2, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 6), mapping="color") +
  theme_bw(base_size = 15) + 
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) + 
  labs(y = "P. multocida- Rel.Abun (All Samples)", x = "Health Status") +
  theme(legend.position="none")

#Date plots
a <- ggplot(data = metadata, aes(x = Date.Collection, y = Mbovis_RelAbunBacteria)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("M. bovis- Rel.Abun (All Samples)") +xlab ("Days relative to end of the study") +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14)) 


b <- ggplot(data = metadata, aes(x = Date.Collection, y = Pm_RelAbunBacteria)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("P. multocida- Rel.Abun (All Samples)") +xlab ("Days relative to end of the study") +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14)) 

ggarrange(a,b, labels = c("A", "B", "C", "D"),
          ncol = 2)

#checking assumptions
# way 1:
plot(m1$full_model)
plot(m2$full_model)
plot(m3$full_model)
plot(m4$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(m1$full_model)) #weird
qqnorm(residuals(m2$full_model)) #weird
qqnorm(residuals(m3$full_model))
qqnorm(residuals(m4$full_model)) #weird

## Only animals that tested positive
str(metadata)
MbYes <- subset(metadata, subset=MbPre_Abs %in% c("1"))
str(MbYes)

HsYes <- subset(metadata, subset=HsPres_Abs %in% c("1"))
str(HsYes)

MhYes <- subset(metadata, subset=MhPres_Abs %in% c("1"))
str(MhYes)

PmYes <- subset(metadata, subset=PmPres_Abs %in% c("1"))
str(MhYes)

#Copy number
y1 <- mixed(Mbovis_copies ~ BRD + Date.Collection + Age + (1|PenCode), data = MbYes,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(y1)
y2 <- mixed(Pm_copies ~ BRD + Date.Collection + Age + (1|PenCode), data = PmYes,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(y2)
y3 <- mixed(Hs_copies ~ BRD + Date.Collection + Age + (1|PenCode), data = HsYes,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(y3)
y4 <- mixed(Mh_copies ~ BRD + Date.Collection + Age + (1|PenCode), data = MhYes,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(y4)

#checking assumptions
# way 1:
plot(y1$full_model)
plot(y2$full_model)
plot(y3$full_model)
plot(y4$full_model)
# this is for testing the normality of the residuals
qqnorm(residuals(y1$full_model)) #weird
qqnorm(residuals(y2$full_model)) #weird
qqnorm(residuals(y3$full_model))
qqnorm(residuals(y4$full_model)) #weird

str(metadata)
y1.1 <- mixed(Mbovis_log ~ BRD + Date.Collection + Age + (1|PenCode), data = MbYes,method = "KR", 
              control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(y1.1)

y1.2 <- mixed(Mbovis_sqrt ~ BRD + Date.Collection + Age + (1|PenCode), data = MbYes,method = "KR", 
              control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(y1.2)

y2.1 <- mixed(Pm_log ~ BRD + Date.Collection + Age + (1|PenCode), data = PmYes,method = "KR", 
              control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(y2.1)

y2.2 <- mixed(Pm_sqrt ~ BRD + Date.Collection + Age + (1|PenCode), data = PmYes,method = "KR", 
              control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(y2.2)


y4.1 <- mixed(Mh_log ~ BRD + Date.Collection + Age + (1|PenCode), data =MhYes,method = "KR", 
              control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(y4.1)

y4.2 <- mixed(Mh_sqrt ~ BRD + Date.Collection + Age + (1|PenCode), data =MhYes,method = "KR", 
              control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(y4.2)


#checking assumptions
# way 1:
plot(y1$full_model)
plot(y1.1$full_model) #better
plot(y1.2$full_model)

plot(y2$full_model)
plot(y2.1$full_model)
plot(y2.2$full_model) #better 

plot(y4$full_model)
plot(y4.1$full_model) #better
plot(y4.2$full_model)  

# this is for testing the normality of the residuals
qqnorm(residuals(y1$full_model)) 
qqnorm(residuals(y1.1$full_model)) #better log
qqnorm(residuals(y1.2$full_model))

qqnorm(residuals(y2$full_model)) 
qqnorm(residuals(y2.1$full_model)) #better log
qqnorm(residuals(y2.2$full_model)) 

qqnorm(residuals(y4$full_model)) 
qqnorm(residuals(y4.1$full_model)) #better log
qqnorm(residuals(y4.2$full_model)) 

##log values are better
#Models to use 
y1.1 
anova(y1.1)
emm_y1.1<- emmeans(y1.1, "BRD")
emm_y1.1

y2.1 
anova(y2.1)

y3
anova(y3)
emm_y3<- emmeans(y3, "BRD")
emm_y3

y4.1 
anova(y4.1)
emm_y4.1<- emmeans(y4.1, "BRD")
emm_y4.1

# interpreting results
afex_plot(y1.1, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 6), mapping="color") +
  theme_bw(base_size = 15) + 
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +  
  labs(y = "M. bovis-(log10)- Positive", x = "Health Status") +
  theme(legend.position="none") 

afex_plot(y3, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 6), mapping="color") +
  theme_bw(base_size = 15) + 
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) + 
  labs(y = "H. somni copies - Positive", x = "Health Status") +
  theme(legend.position="none")

afex_plot(y4.1, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 6), mapping="color") +
  theme_bw(base_size = 15) + 
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +
  labs(y = "M. haemolytica-(log10)-Positive", x = "Health Status") +
  theme(legend.position="none")

c <- ggplot(data = MbYes, aes(x = Date.Collection, y = Mbovis_log)) + 
  geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
  theme_classic() +  ylab("M. bovis-Rel.Abund (log10)-Positive") +xlab ("Days relative to end of the study") +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), axis.title.y = element_text(color="black", size=14, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14)) 

ggarrange(a,c,b, labels = c("A", "B", "C"),
          ncol = 3)

## Relative abundance
#Copy nubmer
str(MbYes)
Y1 <- mixed(Mbovis_RelAbunBacteria ~ BRD + Date.Collection + Age + (1|PenCode), data = MbYes,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(Y1)

str(PmYes)
Y2 <- mixed(Pm_RelAbunBacteria ~ BRD + Date.Collection + Age + (1|PenCode), data = PmYes,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(Y2)

str(HsYes)
Y3 <- mixed(Hs_RelAbunBacteria ~ BRD + Date.Collection + Age + (1|PenCode), data = HsYes,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(Y3)

str(MhYes)
Y4 <- mixed(Mh_RelAbunBacteria ~ BRD + Date.Collection + Age + (1|PenCode), data = MhYes,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(Y4)

#Plots
afex_plot(Y1, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 6), mapping="color") +
  theme_bw(base_size = 15) + 
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) +  
  labs(y = "M. bovis. Rel.Abund-(log10)- Positive", x = "Health Status") +
  theme(legend.position="none") 

afex_plot(Y2, x = "BRD", id = "PenCode", dodge = 0.4, point_arg = list(size = 6), mapping="color") +
  theme_bw(base_size = 15) + 
  theme(legend.position="bottom", panel.grid.major.x = element_blank()) + 
  labs(y = "P. multocida. Rel.Abund-(log10)- Positive", x = "Health Status") +
  theme(legend.position="none")

### Average temperature and relative abundance
## testing temperature
str(metadata)
set_sum_contrasts() 
M1 <- mixed(Mbovis_RelAbunBacteria ~ Ave.temp + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1)

M2 <- mixed(Pm_RelAbunBacteria ~ Ave.temp + (1|PenCode), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)

MbYes <- subset(meta, subset=MbPre_Abs %in% c("1"))
str(MbYes)

M3 <- mixed(Mbovis_RelAbunBacteria ~ Ave.temp + (1|PenCode), data = MbYes,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)

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

