#### NALLEY ET AL. 2022 NUTRIENT THRESHOLD ANALYSIS #### 
### Mortality (Adult Survival, Larval Survival,  Fertilization) ###
rm(list= ls())
library(tidyverse)
library(ggplot2)
library(mixmeta)
library(dosresmeta)
library(splines)
library(effsize)
library(lme4)
library(MuMIn)
library(RColorBrewer)
library(nlme)
library(mgcv)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(ggthemes)
library(mathjaxr)
library(psychmeta)

breaks_log10 <- function(x) {
  low <- floor(log10(min(x)))
  high <- ceiling(log10(max(x)))
  
  10^(seq.int(low, high))
}

mort <- as.data.frame(read.csv("ALL-mortality.csv", header=TRUE))
str(mort) ## everything is factor vs. numeric
mort$SD <- as.numeric(as.character(mort$SD))
mort$SE <- as.numeric(as.character(mort$SE))

# creating new variable that has both Genus and species
mort <- mort %>%
  mutate(Gsp = paste(Genus, Species, sep = "_")) 
head(mort)
mort$RRsig <- as.factor(ifelse(mort$RiskRatio > 1, 1, 0))
mort$lnRR <- log(mort$RiskRatio)
mort$lnSERR <- log(mort$SERiskRatio)
## will give errors if lnSERR = 0, so changing those values to -0.01
mort$lnSERR[mort$lnSERR == -Inf] <- -0.01
mort$type <- "ir"
### adjusting for studies that don't list min value
mort$Nadj <- ifelse(mort$Total.N == 0, 0.1, mort$Total.N)
mort$Padj <- ifelse(mort$PO4 == 0, 0.02, mort$PO4)
### adding log variables
mort$log10N <- log10(mort$Nadj)
mort$log10P <- log10(mort$Padj)
### separating out by age class for later
mort.egg <- as.data.frame(subset(mort, Coral.age.class == "egg"))
mort.larva <- mort[which(mort$Coral.age.class == "larva"),]
mort.adult <- mort[which(mort$Coral.age.class == "adult"),]

## ordering data first
mort.ord <- mort[order(mort$RefID, mort$Control),]
mort.ord.egg <- mort.egg[order(mort.egg$RefID, mort.egg$Control),]
mort.ord.larva <- mort.larva[order(mort.larva$RefID, mort.larva$Control),]
mort.ord.adult <- mort.adult[order(mort.adult$RefID, mort.adult$Control),]

### generating covariance matrix
addS <- lapply(split(mort.ord, mort.ord$RefID), function(x)
  covar.logrr(y=lnRR, v=lnSERR^2, cases=Response, n=N, type=type, data=x))
str(addS)

## egg
addS.egg <- list(addS$NU33a, addS$NU33b, addS$NU33c, addS$NU33d,
                 addS$NU33e, addS$NU33f, addS$NU33g, addS$NU33h,
                 addS$NU33i, addS$NU33j, addS$NU33k, addS$NU33l,
                 addS$NU40a, addS$NU40b, addS$NU40c, addS$NU48a)
names(addS.egg) = c("NU33a", "NU33b", "NU33c", "NU33d", 
                    "NU33e", "NU33f", "NU33g", "NU33h",
                    "NU33i", "NU33j", "NU33k", "NU33l",
                    "NU40a", "NU40b", "NU40c", "NU48a")
## larva
addS.larva <- list(addS$NU02a, addS$NU02b, addS$NU02c, addS$NU11b, addS$NU11d,
                   addS$NU43a, addS$NU43b, addS$NU43c, addS$NU48b, addS$NU48c,
                   addS$NU64a, addS$NU64b, addS$NU73a, addS$NU73b, addS$NU73c,
                   addS$NU73d, addS$NU73e, addS$NU73f)
names(addS.larva) = c("NU02a", "NU02b", "NU02c", "NU11b", "NU11d",
                      "NU43a", "NU43b", "NU43c", "NU48c", "NU48c",
                      "NU64a", "NU64b", "NU73a", "NU73b", "NU73c",
                      "NU73d", "NU73e", "NU73f")

## adult
addS.adult <- list(addS$NU23a, addS$NU23b, addS$NU44, addS$NU47a, addS$NU47b,
                   addS$NU47c, addS$NU63c, addS$NU71)
names(addS.adult) = c("NU23a", "NU23b", "NU44", "NU47a", "NU47b",
                      "NU47c", "NU63c", "NU71")

### removing controls for models
mort2 <- subset(mort, Control=="exp")
mort3 <- subset(mort, SERiskRatio < 5)
mort2$lnRRpos <- as.factor(ifelse(mort2$lnRR > 0, "Increase", "Decrease"))

#
mort2.egg <- as.data.frame(subset(mort2, Coral.age.class == "egg"))
mort2.larva <- as.data.frame(subset(mort2, Coral.age.class == "larva"))
mort2.adult <- as.data.frame(subset(mort2, Coral.age.class == "adult"))

MortTreatSummary <- mort2 %>% 
  group_by(Coral.age.class) %>% 
  summarise(Nmin = min(Nadj),
            Nmax = max(Nadj),
            Pmin = min(Padj),
            Pmax = max(Padj),
            DurMean = mean(Exposure.in.days),
            DurMin = min(Exposure.in.days),
            DurMax = max(Exposure.in.days),
            Exp = n_distinct(RefID),
            Ref = n_distinct(Author.s..Year))

### generating regular effect sizes (SMD) for egg & larva
covar.egg <- by(mort.ord.egg, mort.ord.egg$RefID, function(x) 
  covar.smd(Response, SD, N, "smd", method="hedges", data = x))
covar.egg2 <- list(covar.egg$NU33a, covar.egg$NU33b, covar.egg$NU33c, covar.egg$NU33d,
                  covar.egg$NU33e, covar.egg$NU33f, covar.egg$NU33g, covar.egg$NU33h,
                  covar.egg$NU33i, covar.egg$NU33j, covar.egg$NU33k, covar.egg$NU33l,
                  covar.egg$NU40a, covar.egg$NU40b, covar.egg$NU40c, covar.egg$NU48a)
names(covar.egg2) = c("NU33a", "NU33b", "NU33c", "NU33d", 
                    "NU33e", "NU33f", "NU33g", "NU33h",
                    "NU33i", "NU33j", "NU33k", "NU33l",
                    "NU40a", "NU40b", "NU40c", "NU48a")
mort.ord.egg$smd <- unlist(lapply(covar.egg2, function(x) x$y))
mort.ord.egg$vmd <- unlist(lapply(covar.egg2, function(x) x$v))
str(mort.ord.egg) 
test <- mort.ord.egg[,c("RefID", "Control", "smd")]
mort2b.egg <- subset(mort.ord.egg, Control == "exp")
mort2b.egg$smdpos <- as.factor(ifelse(mort2b.egg$smd > 0, "Increase", "Decrease"))
newlist.egg <- list(NA)
for (i in seq(1,length(covar.egg2))) {
  newlist.egg[i] <- list(covar.egg2[[i]]$S)
}
##
# LINEAR FIXED AND RANDOM, NESTED EFFECTS ACCOUNTING FOR WITHIN-COMPARISON & WITHIN-STUDY CORRELATIONS
#Code for covariance matrix to include for addSlist for Ref/Comparison
newlist2.egg <- list(NA) #collects the list of covariance matrices for the block diag matrix by reference for the nested hierarchial ref/comparison
templist.egg <- list(NA) #holds temp list of covariance matrices associated with each reference
templist2.egg <-list(NA)
reflista.egg <- mort.ord.egg %>% distinct(Author.s..Year,RefID)
reflist.egg <- reflista.egg[,2]
for (i in seq(1,length(unique(reflista.egg$Author.s..Year))))  {
  #pull the elements from covar_DS_photo that are all from the same reference [i]
  templist.egg[i] <-list(covar.egg2[reflist.egg==unique(reflista.egg$Author.s..Year)[i]])
  for (j in seq(1,length(templist.egg[[i]]))) {
    #for each comparison in the reference, pull out the covar matrices (element $S) and put in into templist2
    templist2.egg[j] <- list(templist.egg[[i]][[j]]$S)
  }
  #turn list of covars from all comparison in one reference into block diag matrix
  newlist2.egg[i] <- list(bdiagMat(templist2.egg))
  templist2.egg <- list(NA)
}


### generating regular effect sizes (SMD) for egg & larva
unique(mort.ord.larva$RefID)
covar.larva <- by(mort.ord.larva, mort.ord.larva$RefID, function(x) 
  covar.smd(Response, SD, N, "smd", method="hedges", data = x))
covar.larva2 <- list(covar.larva$NU02a, covar.larva$NU02b, covar.larva$NU02c, covar.larva$NU11b, covar.larva$NU11d,
                     covar.larva$NU43a, covar.larva$NU43b, covar.larva$NU43c, covar.larva$NU48b, covar.larva$NU48c,
                     covar.larva$NU64a, covar.larva$NU64b, covar.larva$NU73a, covar.larva$NU73b, covar.larva$NU73c,
                     covar.larva$NU73d, covar.larva$NU73e, covar.larva$NU73f)
names(covar.larva2) = c("NU02a", "NU02b", "NU02c", "NU11b", "NU11d",
                      "NU43a", "NU43b", "NU43c", "NU48b", "NU48c",
                      "NU64a", "NU64b", "NU73a", "NU73b", "NU73c",
                      "NU73d", "NU73e", "NU73f")
mort.ord.larva$smd <- unlist(lapply(covar.larva2, function(x) x$y))
mort.ord.larva$vmd <- unlist(lapply(covar.larva2, function(x) x$v))
str(mort.ord.larva) 
test <- mort.ord.larva[,c("RefID", "Control", "smd")]
mort2b.larva <- subset(mort.ord.larva, Control == "exp")
mort2b.larva$smdpos <- as.factor(ifelse(mort2b.larva$smd > 0, "Increase", "Decrease"))
newlist.larva <- list(NA)
for (i in seq(1,length(covar.larva2))) {
  newlist.larva[i] <- list(covar.larva2[[i]]$S)
}
###
# LINEAR FIXED AND RANDOM, NESTED EFFECTS ACCOUNTING FOR WITHIN-COMPARISON & WITHIN-STUDY CORRELATIONS
#Code for covariance matrix to include for addSlist for Ref/Comparison
newlist2 <- list(NA) #collects the list of covariance matrices for the block diag matrix by reference for the nested hierarchial ref/comparison
templist <- list(NA) #holds temp list of covariance matrices associated with each reference
templist2 <-list(NA)
reflista <- mort.ord.larva %>% distinct(Author.s..Year,RefID)
reflist <- reflista[,2]
for (i in seq(1,length(unique(reflista$Author.s..Year))))  {
  #pull the elements from covar_DS_photo that are all from the same reference [i]
  templist[i] <-list(covar.larva2[reflist==unique(reflista$Author.s..Year)[i]])
  for (j in seq(1,length(templist[[i]]))) {
    #for each comparison in the reference, pull out the covar matrices (element $S) and put in into templist2
    templist2[j] <- list(templist[[i]][[j]]$S)
  }
  #turn list of covars from all comparison in one reference into block diag matrix
  newlist2[i] <- list(bdiagMat(templist2))
  templist2 <- list(NA)
}



## dotplot
## adult
ggplot(mort2.adult, 
       aes(x = Nadj, y = Padj, 
           color = lnRRpos,
           size = abs(lnRR))) + 
  ### surface water from station aloha 2019 annual report (surface ocean)
  annotate("point", x = 0.03, y = 0.03, colour = "gray", size = 4, shape = 8) +
  ## ambient reef from Nakajima et al. 2015 (Malaysia)
  annotate("point", x = 0.75, y = 0.1, colour = "gray", size = 4, shape = 8) +
  ## ambient reef from Silbiger et al. 2018 (Kaneohe Bay)
  annotate("point", x = .15, y = .15, colour = "gray", size = 4, shape = 8) + 
  ## ambient (but high) from AIMS reef in Fabricius et al. 2013 (GBR)
  annotate("point", x = .44, y = .38, colour = "gray", size = 4, shape = 8) + 
  ## high reef from Silbiger et al. 2018 (throughout Pacific)
  annotate("point", x = 7.6, y = 2.6, colour = "gray", size = 4, shape = 8) + 
  ## highest measured at SGD seep in Lubarsky et al. 2018 (Maunalua Bay)
  annotate("point", x = 32.39, y = 0.89, colour = "gray", size = 4, shape = 8) + 
  scale_x_log10(breaks = breaks_log10) +  
  scale_y_log10(breaks = breaks_log10) + annotation_logticks() +
  geom_jitter(alpha = .5, position = position_jitter(width = .07, height = 0.03)) +
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  scale_size(range = c(3,15)) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = expression("DIP ("*mu*"M)"),
       shape = "Age Class",
       color = "% Survival",
       size = "Treatment\nEffect Size") + 
  guides(color = guide_legend(order = 2, override.aes = list(size=4)),
         shape = guide_legend(order = 3, override.aes = list(size=4)),
         size = guide_legend(order = 1)) 

## larva
ggplot(mort2b.larva, 
       aes(x = Nadj, y = Padj, 
           color = smdpos,
           size = abs(smd))) + 
  ## surface water from station aloha 2019 annual report (surface ocean)
  annotate("point", x = 0.03, y = 0.03, colour = "gray", size = 4, shape = 8) +
  ## ambient reef from Nakajima et al. 2015 (Malaysia)
  annotate("point", x = 0.75, y = 0.1, colour = "gray", size = 4, shape = 8) +
  ## ambient reef from Silbiger et al. 2018 (Kaneohe Bay)
  annotate("point", x = .15, y = .15, colour = "gray", size = 4, shape = 8) + 
  ## ambient (but high) from AIMS reef in Fabricius et al. 2013 (GBR)
  annotate("point", x = .44, y = .38, colour = "gray", size = 4, shape = 8) + 
  ## high reef from Silbiger et al. 2018 (throughout Pacific)
  annotate("point", x = 7.6, y = 2.6, colour = "gray", size = 4, shape = 8) + 
  ## highest measured at SGD seep in Lubarsky et al. 2018 (Maunalua Bay)
  annotate("point", x = 32.39, y = 0.89, colour = "gray", size = 4, shape = 8) + 
  scale_x_log10(breaks = breaks_log10) +  
  scale_y_log10(breaks = breaks_log10) + annotation_logticks() +
  geom_jitter(alpha = .5, position = position_jitter(width = .07, height = 0.03)) +
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  scale_size(range = c(3,15)) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = expression("DIP ("*mu*"M)"),
       shape = "Age Class",
       color = "% Survival",
       size = "Treatment\nEffect Size") + 
  guides(color = guide_legend(order = 2, override.aes = list(size=4)),
         shape = guide_legend(order = 3, override.aes = list(size=4)),
         size = guide_legend(order = 1)) 

## egg
ggplot(mort2b.egg, 
       aes(x = Nadj, y = Padj, 
           color = smdpos,
           size = abs(smd))) + 
  ## surface water from station aloha 2019 annual report (surface ocean)
  annotate("point", x = 0.03, y = 0.03, colour = "gray", size = 4, shape = 8) +
  ## ambient reef from Nakajima et al. 2015 (Malaysia)
  annotate("point", x = 0.75, y = 0.1, colour = "gray", size = 4, shape = 8) +
  ## ambient reef from Silbiger et al. 2018 (Kaneohe Bay)
  annotate("point", x = .15, y = .15, colour = "gray", size = 4, shape = 8) + 
  ## ambient (but high) from AIMS reef in Fabricius et al. 2013 (GBR)
  annotate("point", x = .44, y = .38, colour = "gray", size = 4, shape = 8) + 
  ## high reef from Silbiger et al. 2018 (throughout Pacific)
  annotate("point", x = 7.6, y = 2.6, colour = "gray", size = 4, shape = 8) + 
  ## highest measured at SGD seep in Lubarsky et al. 2018 (Maunalua Bay)
  annotate("point", x = 32.39, y = 0.89, colour = "gray", size = 4, shape = 8) + 
  scale_x_log10(breaks = breaks_log10) +  
  scale_y_log10(breaks = breaks_log10) + annotation_logticks() +
  geom_jitter(alpha = .5, position = position_jitter(width = .07, height = 0.03)) +
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  scale_size(range = c(3,15)) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = expression("DIP ("*mu*"M)"),
       shape = "Age Class",
       color = "Fertilization\nSuccess",
       size = "Treatment\nEffect Size") + 
  guides(color = guide_legend(order = 2, override.aes = list(size=4)),
         shape = guide_legend(order = 3, override.aes = list(size=4)),
         size = guide_legend(order = 1)) 

### mort vs. adult
MortN.adult <- ggplot(mort.adult, 
               aes(x = Nadj, y = Response,
                   color = Control)) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  geom_smooth(method = "loess", color = "darkgrey", fill = "lightgrey") +
  geom_point(alpha = 0.7) + theme_bw() + scale_y_continuous(limits = c(-10,120),
                                                            breaks = c(0, 20, 40, 60, 80,100))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.15, 0.95),
        legend.background = element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=16)) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "% Adult Survival",
       color =  element_blank())  +
  scale_color_brewer(palette="Set2") 
MortN.adult
MortP.adult <- ggplot(mort.adult, 
               aes(x = Padj, y = Response,
                   color = Control)) + 
  scale_y_continuous(limits = c(-10,120),
                     breaks = c(0, 20, 40, 60, 80,100)) + scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  geom_smooth(method = "loess", color = "darkgrey", fill = "lightgrey") +
  geom_point(alpha = 0.7) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        text = element_text(size=20),
        legend.position = "none",
        legend.background = element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "% Survival",
       color =  expression("DIN ("*mu*"M)"))  +
  scale_color_brewer(palette="Set2")
MortP.adult
grid.arrange(MortN.adult,MortP.adult,ncol=2)

### mort vs. larva
MortN.larva <- ggplot(mort.larva, 
                      aes(x = Nadj, y = Response,
                          color = Control)) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  geom_smooth(method = "loess", color = "darkgrey", fill = "lightgrey") +
  geom_point(alpha = 0.7) + theme_bw() + ylim(0,100) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.15, 0.85),
        legend.background = element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=12)) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "% Larval Survival",
       color =  element_blank())   +
  scale_color_brewer(palette="Set2")
MortN.larva

MortP.larva <- ggplot(mort.larva, 
                      aes(x = Padj, y = Response,
                          color = Control)) + 
  ylim(0,100) + scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  geom_smooth(method = "loess", color = "darkgrey", fill = "lightgrey") +
  geom_point(alpha = 0.7) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        text = element_text(size=20),
        legend.position = "none",
        legend.background = element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "% Survival",
       color =  expression("Total N "*mu*"M"))   +
  scale_color_brewer(palette="Set2")
MortP.larva
grid.arrange(MortN.larva,MortP.larva,ncol=2)


### mort vs. egg
MortN.egg <- ggplot(mort.egg, 
                      aes(x = Nadj, y = Response,
                          color = Control)) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  geom_smooth(method = "loess", color = "darkgrey", fill = "lightgrey") +
  geom_point(alpha = 0.7) + theme_bw() + ylim(0,100) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.15, 0.2),
        legend.background = element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=12)) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "% Fertilization",
       color = element_blank()) +
  scale_color_brewer(palette="Set2")
MortN.egg

MortP.egg <- ggplot(mort.egg, 
                      aes(x = Padj, y = Response,
                          color = Control)) + 
  ylim(0,100) + scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  geom_smooth(method = "loess", color = "darkgrey", fill = "lightgrey") +
  geom_point(alpha = 0.7) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        text = element_text(size=20),
        legend.position = "none",
        legend.background = element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "% Survival",
       color =  element_blank()) +
  scale_color_brewer(palette="Set2")
MortP.egg
grid.arrange(MortN.egg,MortP.egg,ncol=2)

###### EGG #### 
## fert 
plot(lnRR ~ log10N, data = mort2.egg)
plot(smd ~ log10N, data = mort2b.egg)
plot(lnRR ~ log10P, data = mort2.egg)
plot(smd ~ log10P, data = mort2b.egg)

egg1 <- mixmeta(smd ~ 0 + log10N + log10P,
                    random = ~ 1 | RefID,
                    data=mort2b.egg, method="ml",
                control=list(addSlist=newlist.egg))
summary(egg1)
egg1b <- mixmeta(lnRR ~ 0 + log10N + log10P,
                random = ~ 1 | RefID,
                data=mort2.egg, method="ml",
                control=list(addSlist=addS.egg))
summary(egg1b)
egg2 <- mixmeta(smd ~ log10N + log10P,
                random = ~ 1 | RefID,
                data=mort2b.egg, method="ml",
                control=list(addSlist=newlist.egg))
summary(egg2)
egg2b <- mixmeta(lnRR ~ log10N + log10P,
                random = ~ 1 | RefID,
                data=mort2.egg, method="ml",
                control=list(addSlist=addS.egg))
summary(egg2b)
egg2c <- mixmeta(smd ~ log10N + log10P,
                random = ~ 1 |Author.s..Year/RefID,
                data=mort2b.egg, method="ml",
                control=list(addSlist=newlist2.egg))
summary(egg2c)
egg2d <- mixmeta(lnRR ~ log10N + log10P ,
                 random =  list(~1|Gsp, ~ 1 | RefID),
                 data=mort2.egg, method="ml",
                 S = lnSERR^2)
summary(egg2d)

egg3 <- mixmeta(smd ~ 0 + log10N + log10P + Exposure.in.days,
                random = ~ 1 | RefID,
                data=mort2b.egg, method="ml",
                control=list(addSlist=newlist.egg))
summary(egg3)
egg3b <- mixmeta(lnRR ~ 0 + log10N + log10P + Exposure.in.days,
                random = ~ 1 | RefID,
                data=mort2.egg, method="ml",
                control=list(addSlist=addS.egg))
summary(egg3b)
egg4 <- mixmeta(smd ~ Nadj + Padj,
                random = ~ 1 | RefID,
                data=mort2b.egg, method="ml",
                control=list(addSlist=newlist.egg))
summary(egg4)
egg4b <- mixmeta(lnRR ~ 0 + Nadj + Padj,
                 random = ~ 1 | RefID,
                 data=mort2.egg, method="ml",
                 control=list(addSlist=addS.egg))
summary(egg4b)
egg4c <- mixmeta(smd ~ Nadj + Padj,
                 random = ~ 1 | Author.s..Year/RefID,
                 data=mort2b.egg, method="ml",
                 control=list(addSlist=newlist2.egg))
summary(egg4c)
egg5 <- mixmeta(smd ~ 0 + Nadj + Padj,
                random = ~ 1 | RefID,
                data=mort2b.egg, method="ml",
                control=list(addSlist=newlist.egg))
summary(egg5)
egg5b <- mixmeta(lnRR ~ 0 + Nadj + Padj,
                 random = ~ 1 | RefID,
                 data=mort2.egg, method="ml",
                 control=list(addSlist=addS.egg))
summary(egg5b)
egg6 <- mixmeta(smd ~ Nadj + Padj + Exposure.in.days,
                random = ~ 1 | RefID,
                data=mort2b.egg, method="ml",
                control=list(addSlist=newlist.egg))
summary(egg6)
egg6b <- mixmeta(lnRR ~ Nadj + Padj + Exposure.in.days,
                 random = ~ 1 | RefID,
                 data=mort2.egg, method="ml",
                 control=list(addSlist=addS.egg))
summary(egg6b)
mort2.egg$N.exp <- mort2.egg$Nadj/mort2.egg$Exposure.in.days
mort2.egg$P.exp <- mort2.egg$Padj/mort2.egg$Exposure.in.days
egg6c <- mixmeta(lnRR ~ N.exp + P.exp,
                 random = ~ 1 | RefID,
                 data=mort2.egg, method="ml",
                 control=list(addSlist=addS.egg))
summary(egg6c)

AIC(egg1, egg1b, egg2, egg2b,egg2d, egg3, egg3b, egg4, egg4b, egg4c, egg5, egg5b, egg6, egg6b, egg6c)
BIC(egg1, egg1b, egg2, egg2b, egg2d,egg3, egg3b, egg4, egg4b, egg5, egg5b, egg6, egg6b)
summary(egg2b)
summary(egg6b)

## looking at best fit with smd not lnRR
resid <- resid(egg4c)
fitted <- fitted(egg4c)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)

### looking at nonlinear
## generating 1st & 3rd Q for nonlinear models
summary(mort2.egg$Nadj)
summary(mort2.egg$Padj)
summary(mort2.egg$log10N)
summary(mort2.egg$log10P)

## simple random
mod7 <- mixmeta(smd ~ ns(log10N, knots=c(0.02788,1.32804)) + 
                   ns(log10P, knots=c(-0.69897,0.63055)), 
                 data=mort2b.egg,
                 random= ~ 1| RefID, 
                 method="ml", 
                 control=list(addSlist=newlist.egg))
summary(mod7) #
mod7b <- mixmeta(lnRR ~ ns(log10N, knots=c(0.02788,1.32804)) + 
                  ns(log10P, knots=c(-0.69897,0.63055)), 
                data=mort2.egg,
                random= ~ 1| RefID, 
                method="ml", 
                control=list(addSlist=addS.egg))
summary(mod7b)
mod7c <- mixmeta(lnRR ~ ns(log10N, knots=c(0.02788,1.32804)) + 
                   ns(log10P, knots=c(-0.69897,0.63055)), 
                 data=mort2.egg,
                 random= list(~1|Gsp),
                 S=lnSERR^2)
summary(mod7c)
mod8 <- mixmeta(smd ~ ns(Nadj, knots=c(1.075,21.413)) + 
                  ns(Padj, knots=c(0.20,4.30)), 
                data=mort2b.egg,
                random= ~ 1| RefID, 
                method="ml", 
                control=list(addSlist=newlist.egg))
summary(mod8) #
mod8b <- mixmeta(lnRR ~ ns(Nadj, knots=c(1.075,21.413)) + 
                   ns(Padj, knots=c(0.20,4.30)), 
                 data=mort2.egg,
                 random= ~ 1| RefID, 
                 method="ml", 
                 control=list(addSlist=addS.egg))
summary(mod8b)
## with quartiles
mod9 <- mixmeta(smd ~ ns(log10N, knots=c(0.02788,0.59864,1.32804)) + 
                  ns(log10P, knots=c(-0.69897,-0.15490,0.63055)), 
                data=mort2b.egg,
                random= ~ 1| RefID, 
                method="ml", 
                control=list(addSlist=newlist.egg))
summary(mod9) #
mod9b <- mixmeta(lnRR ~ ns(log10N, knots=c(0.02788,0.59864,1.32804)) + 
                   ns(log10P, knots=c(-0.69897,-0.15490,0.63055)), 
                 data=mort2b.egg,
                 random= ~ 1| RefID, 
                 method="ml", 
                 control=list(addSlist=addS.egg))
summary(mod9b)
mod10 <- mixmeta(smd ~ ns(Nadj, knots=c(1.075,4.075,21.413)) + 
                   ns(Padj, knots=c(0.20,0.70,4.30)), 
                 data=mort2b.egg,
                 random= ~ 1| RefID, 
                 method="ml", 
                control=list(addSlist=newlist.egg))
summary(mod10) #
mod10b <- mixmeta(lnRR ~ ns(Nadj, knots=c(1.075,4.075,21.413)) + 
                   ns(Padj, knots=c(0.20,0.70,4.30)), 
                 data=mort2b.egg,
                 random= ~ 1| RefID, 
                 method="ml", 
                 control=list(addSlist=addS.egg))
summary(mod10b)

AIC(egg2b, egg2d, egg6b, egg4, mod7, mod7b, mod7c, mod8, mod8b, mod9, mod9b, mod10, mod10b)
BIC(egg2b, egg6b, egg4, mod7, mod7b, mod8, mod8b, mod9, mod9b, mod10, mod10b)
##
summary(egg2d)
resid <- resid(egg2d)
fitted <- fitted(egg2d)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)

### plotting with predictions
newdata.predictN.egg <- data.frame(log10N = mort2b.egg$log10N,
                               log10P = median(mort2b.egg$log10P))
newdata.predictP.egg <- data.frame(log10N = median(mort2b.egg$log10N),
                               log10P = mort2b.egg$log10P)
predN.egg <- predict(egg2d, newdata=newdata.predictN.egg, ci=TRUE)
predP.egg <- predict(egg2d, newdata=newdata.predictP.egg, ci=TRUE)
testCIN.egg <- cbind(mort2b.egg, predN.egg)
testCIP.egg <- cbind (mort2b.egg, predP.egg)
min_conc_N.egg <- testCIN.egg %>% 
  mutate(overlap0 = 0 >= ci.lb & 0 <= ci.ub) %>% 
  filter(overlap0==FALSE) %>% 
  summarize(min_N = min(Nadj),
            max_N = max(Nadj),
            min_P = min(Padj),
            max_P = max(Padj))

## acropora longicyathus was much lower than the others
## and now plotting for nitrogen as a predictor
ESNegg <- ggplot(testCIN.egg, 
                 aes(x = Nadj,
                     y = lnRR, color = Padj,
                     ymin = lnRR-lnSERR,
                     ymax = lnRR+lnSERR)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Log Risk Ratio (+/- s.e.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.85, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) +  
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIP ("*mu*"M)")) 
ESNegg

ESPegg <- ggplot(testCIP.egg, 
                 aes(x = Padj,
                     y = lnRR, color = Nadj,
                     ymin = lnRR-lnSERR,
                     ymax = lnRR+lnSERR)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)")) +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = c(0.15, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIN ("*mu*"M)")) +
  geom_ribbon(aes(x = Padj, y = lnRR,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Padj, y = fit), inherit.aes=FALSE) 
ESPegg2 <- ggplot(testCIP.egg, 
                 aes(x = Padj,
                     y = lnRR, color = Gsp,
                     ymin = lnRR-lnSERR,
                     ymax = lnRR+lnSERR)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Log Risk Ratio (+/- s.e.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
       legend.position = c(0.15, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("Species")) 
ESPegg2
grid.arrange(ESNegg,ESPegg,ncol=2)
summary(egg2d)

### quartile plots
eggplot3 <- mort2b.egg%>%mutate(NadjQuantBins = cut(Nadj, 
                                                        breaks = unique(quantile(Nadj,probs=seq.int(0,1, by=1/4))), 
                                                        include.lowest=TRUE))
eggplot3 <- eggplot3%>%mutate(PadjQuantBins = cut(Padj, 
                                                      breaks = unique(quantile(Padj,probs=seq.int(0,1, by=1/4))), 
                                                      include.lowest=TRUE))
N <- ggplot(eggplot3, aes(x = NadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Proportion of Reported Change in Fertilization",
       fill = "") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(size = 15)) +
  scale_x_discrete(labels = c("[0.65,1.07]" = "0.65-1.07",
                              "(1.07,4.08]" = "1.07-4.08",
                              "(4.08,21.4]" = "4.08-21.4",
                              "(21.4,202]" = "21.4-202"))
N
P <- ggplot(eggplot3, aes(x = PadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "",
       fill = "") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 15))  +
  scale_x_discrete(labels = c("[0.08,0.2]" = "0.08-0.2",
                              "(0.2,0.7]" = "0.2-0.7",
                              "(0.7,4.3]" = "0.7-4.3",
                              "(4.3,101]" = "4.3-101"))
P
grid.arrange(N,P,ncol=2)

###### LARVA ####
## larva 
larva1 <- mixmeta(smd ~ log10N + log10P,
                random = ~ 1 |RefID,
                data=mort2b.larva, method="ml",
                control=list(addSlist=newlist.larva))
summary(larva1)
larva1b <- mixmeta(lnRR ~ log10N + log10P,
                 random = ~ 1 | RefID,
                 data=mort2b.larva, method="ml",
                 control=list(addSlist=addS.larva))
summary(larva1b)
larva1c <- mixmeta(smd ~ log10N + log10P,
                  random = ~ 1 |Author.s..Year/RefID,
                  data=mort2b.larva, method="ml",
                  control=list(addSlist=newlist2))
summary(larva1c)
AIC(larvax, larva1, larva1c)

larva2 <- mixmeta(smd ~ 0+ log10N + log10P,
                random = ~ 1 | RefID,
                data=mort2b.larva, method="ml",
                control=list(addSlist=newlist.larva))
summary(larva2)
larva2b <- mixmeta(lnRR ~ 0+ log10N + log10P,
                 random = ~ 1 | RefID,
                 data=mort2b.larva, method="ml",
                 control=list(addSlist=addS.larva))
summary(larva2b)
larva3 <- mixmeta(smd ~ 0 + log10N + log10P + Exposure.in.days,
                random = ~ 1 | RefID,
                data=mort2b.larva, method="ml",
                control=list(addSlist=newlist.larva))
summary(larva3)
larva3b <- mixmeta(lnRR ~ 0 + log10N + log10P + Exposure.in.days,
                 random = ~ 1 | RefID,
                 data=mort2b.larva, method="ml",
                 control=list(addSlist=addS.larva))
summary(larva3b)
larva4 <- mixmeta(smd ~ Nadj + Padj,
                random = ~ 1 | RefID,
                data=mort2b.larva, method="ml",
                control=list(addSlist=newlist.larva))
summary(larva4)
larva4b <- mixmeta(lnRR ~ Nadj + Padj,
                 random = ~ 1 | RefID,
                 data=mort2b.larva, method="ml",
                 control=list(addSlist=addS.larva))

summary(larva4b)
larva4c <- mixmeta(smd ~ Nadj + Padj,
                   random = ~ 1 |Author.s..Year/RefID,
                   data=mort2b.larva, method="ml",
                   control=list(addSlist=newlist2))
summary(larva4c)
larva4d <- mixmeta(smd ~ Nadj + Padj,
                   random = ~ 1 |Author.s..Year/RefID,
                   data=mort2b.larva, method="ml",
                   S = lnSERR^2)
summary(larva4d)

larva4e <- mixmeta(lnRR ~ log10N + log10P,
                   random =  list(~1|Gsp, ~ 1 | RefID),
                   data=mort2.larva, method="ml",
                   S = lnSERR^2)
summary(larva4e)

larva5 <- mixmeta(smd ~ 0 + Nadj + Padj,
                random = ~ 1 | RefID,
                data=mort2b.larva, method="ml",
                control=list(addSlist=newlist.larva))
summary(larva5)
larva5b <- mixmeta(lnRR ~ 0 + Nadj + Padj,
                 random = ~ 1 | RefID,
                 data=mort2b.larva, method="ml",
                 control=list(addSlist=addS.larva))
summary(larva5b)

larva6 <- mixmeta(smd ~ Nadj + Padj + Exposure.in.days,
                random = ~ 1 | RefID,
                data=mort2b.larva, method="ml",
                control=list(addSlist=newlist.larva))
summary(larva6)
larva6b <- mixmeta(lnRR ~ Nadj + Padj + Exposure.in.days,
                 random = ~ 1 | RefID,
                 data=mort2b.larva, method="ml",
                 control=list(addSlist=addS.larva))
summary(larva6b)
AIC(larva1, larva1b,larva2, larva2b, larva3, larva3b, larva4, larva4b, larva5, larva5b,larva6, larva6b)
BIC(larva1, larva1b,larva2, larva2b, larva3, larva3b, larva4, larva4b, larva5, larva5b,larva6, larva6b)
summary(larva4b)
summary(larva5b)
resid <- resid(larva5b)
fitted <- fitted(larva5b)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)
## one outlier

### looking at nonlinear
## generating 1st & 3rd Q for nonlinear models
summary(mort2b.larva$Nadj)
summary(mort2b.larva$Padj)
summary(mort2b.larva$log10N)
summary(mort2b.larva$log10P)
plot(lnRR ~ log10N + log10P, data = mort2b.larva)
## simple random
mod7 <- mixmeta(smd ~ ns(log10N, knots=c(1.012,1.718)) + 
                  ns(log10P, knots=c(-1.2812,-0.4292)), 
                data=mort2b.larva,
                random= ~ 1| RefID, 
                method="ml", 
                control=list(addSlist=newlist.larva))
summary(mod7) #
mod7b <- mixmeta(lnRR ~ ns(log10N, knots=c(1.012,1.718)) + 
                   ns(log10P, knots=c(-1.2812,-0.4292)), 
                 data=mort2b.larva,
                 random= ~ 1| RefID, 
                 method="ml", 
                 control=list(addSlist=addS.larva))
summary(mod7b)
mod8 <- mixmeta(smd ~ ns(Nadj, knots=c(10.29,52.25)) + 
                  ns(Padj, knots=c(0.0525,0.3723)), 
                data=mort2b.larva,
                random= ~ 1| RefID, 
                method="ml", 
                control=list(addSlist=newlist.larva))
summary(mod8) #
mod8b <- mixmeta(lnRR ~  ns(Nadj, knots=c(10.29,52.25)) + 
                   ns(Padj, knots=c(0.0525,0.3723)), 
                 data=mort2b.larva,
                 random= ~ 1| RefID, 
                 method="ml", 
                 control=list(addSlist=addS.larva))
summary(mod8b)
## with quartiles
mod9 <- mixmeta(smd ~ ns(log10N, knots=c(1.012,1.277,1.718)) + 
                  ns(log10P, knots=c(-1.2812,-1.0969,0.4292)), 
                data=mort2b.larva,
                random= ~ 1| RefID, 
                method="ml", 
                control=list(addSlist=newlist.larva))
summary(mod9) #
mod9b <- mixmeta(lnRR ~ ns(log10N, knots=c(1.012,1.277,1.718)) + 
                   ns(log10P, knots=c(-1.2812,-1.0969,0.4292)), 
                 data=mort2b.larva,
                 random= ~ 1| RefID, 
                 method="ml", 
                 control=list(addSlist=addS.larva))
summary(mod9b)
mod10 <- mixmeta(smd ~ ns(Nadj, knots=c(10.29,18.95,52.25)) + 
                   ns(Padj, knots=c(0.0525,0.0800,0.3723)), 
                 data=mort2b.larva,
                 random= ~ 1| RefID, 
                 method="ml", 
                 control=list(addSlist=newlist.larva))
summary(mod10) #
mod10b <- mixmeta(lnRR ~ ns(Nadj, knots=c(10.29,18.95,52.25)) + 
                    ns(Padj, knots=c(0.0525,0.0800,0.3723)), 
                  data=mort2b.larva,
                  random= ~ 1| RefID, 
                  method="ml", 
                  control=list(addSlist=addS.larva))
summary(mod10b)
summary(larva5b)
AIC(larva4, larva4e, larva5b, mod7, mod7b, mod8, mod8b, mod9,mod9b, mod10, mod10b)
BIC(larva4, larva4e, larva5b, mod7, mod7b, mod8, mod8b, mod9,mod9b, mod10, mod10b)

summary(larva4e)
resid <- resid(larva4e)
fitted <- fitted(larva4e)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)
## one outlier

summary(mod7b)
resid <- resid(mod7b)
fitted <- fitted(mod7b)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)
## one outlier

### plotting with predictions
newdata.predictN.larva <- data.frame(log10N = mort2b.larva$log10N,
                                   log10P = median(mort2b.larva$log10P))
newdata.predictP.larva <- data.frame(log10N = median(mort2b.larva$log10N),
                                   log10P = mort2b.larva$log10P)
predN.larva <- predict(larva4e, newdata=newdata.predictN.larva, ci=TRUE)
predP.larva <- predict(larva4e, newdata=newdata.predictP.larva, ci=TRUE)
testCIN.larva <- cbind(mort2b.larva, predN.larva)
testCIP.larva <- cbind (mort2b.larva, predP.larva)
min_conc_N.larva <- testCIN.larva %>% 
  mutate(overlap0 = 0 >= ci.lb & 0 <= ci.ub) %>% 
  filter(overlap0==FALSE) %>% 
  summarize(min_N = min(Nadj),
            max_N = max(Nadj),
            min_P = min(Padj),
            max_P = max(Padj))

## and now plotting for nitrogen as a predictor
ESNlarva <- ggplot(testCIN.larva, 
                 aes(x = Nadj,
                     y = lnRR, color = Padj,
                     ymin = lnRR-lnSERR,
                     ymax = lnRR+lnSERR)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Log Risk Ratio (+/- s.e.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.15, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) +  
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = Nadj, y = lnRR,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) 
ESNlarva
ESPlarva <- ggplot(testCIP.larva, 
                 aes(x = Padj,
                     y = lnRR, color = Nadj,
                     ymin = lnRR-lnSERR,
                     ymax = lnRR+lnSERR)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)")) +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = c(0.85, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIN ("*mu*"M)"))  
grid.arrange(ESNlarva,ESPlarva,ncol=2)
summary(larva4e)

### quartile plots
larvaplot3 <- mort2b.larva%>%mutate(NadjQuantBins = cut(Nadj, 
                                                breaks = unique(quantile(Nadj,probs=seq.int(0,1, by=1/4))), 
                                                include.lowest=TRUE))
larvaplot3 <- larvaplot3%>%mutate(PadjQuantBins = cut(Padj, 
                                                breaks = unique(quantile(Padj,probs=seq.int(0,1, by=1/4))), 
                                                include.lowest=TRUE))
N <- ggplot(larvaplot3, aes(x = NadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Proportion of Reported Change in Larval Survival",
       fill = "") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(size = 15)) +
  scale_x_discrete(labels = c("[0.06,10.3]" = "0.06-10.3",
                              "(10.3,18.9]" = "10.3-18.9",
                              "(18.9,52.2]" = "18.9-52.2",
                              "(52.2,202]" = "52.2-202"))
P <- ggplot(larvaplot3, aes(x = PadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "",
       fill = "") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 15))  +
  scale_x_discrete(labels = c("[0.02,0.0525]" = "0.02-0.05",
                              "(0.0525,0.08]" = "0.05-0.08",
                              "(0.08,0.372]" = "0.08-0.37",
                              "(0.372,100]" = "0.37-100"))
grid.arrange(N,P,ncol=2)



##### adult ####
adult1 <- mixmeta(lnRR ~  log10N + log10P,
                  random = ~ 1 | RefID,
                  data=mort2.adult, method="ml",
                  control=list(addSlist=addS.adult))
summary(adult1)
adult1b <- mixmeta(lnRR ~ Nadj + Padj,
                  random = ~ 1 | RefID,
                  data=mort2.adult, method="ml",
                  control=list(addSlist=addS.adult))
summary(adult1b)
AIC(adult1b, adult1)
BIC(adult1b, adult1)
## adult1

adult1c <- mixmeta(lnRR ~ Nadj + Padj + Nadj:Padj,
                   random = ~ 1 | RefID,
                   data=mort2.adult, method="ml",
                   control=list(addSlist=addS.adult))
summary(adult1c)

adult1d <- mixmeta(lnRR ~ Nadj + Padj + Exposure.in.days,
                   random = ~ 1 | RefID,
                   data=mort2.adult, method="ml",
                   control=list(addSlist=addS.adult))
summary(adult1d)

adult1e <- mixmeta(lnRR ~ Nadj + Padj + Gsp,
                   random = ~ 1 | RefID,
                   data=mort2.adult, method="ml",
                   control=list(addSlist=addS.adult))
summary(adult1e)

adult3 <- mixmeta(lnRR ~ 0 + log10N + log10P + Exposure.in.days,
                  random = ~ 1 | RefID,
                  data=mort2.adult, method="ml",
                  control=list(addSlist=addS.adult))
summary(adult3)

adult4 <- mixmeta(lnRR ~ log10N + log10P, S = lnSERR^2, data = mort2.adult,
                  random = list(~1 | Gsp, ~1 | RefID), method = "ml")
summary(adult4)

adult5 <- mixmeta(lnRR ~ Nadj + Padj,
                  random =  list(~1|Gsp, ~ 1 | RefID),
                  data=mort2.adult, method="ml",
                  S = lnSERR^2)
summary(adult5)
adult5b <- mixmeta(lnRR ~ Nadj + Padj,
                  random =  ~1|Gsp/RefID, 
                  data=mort2.adult, method="ml",
                  S = lnSERR^2)
summary(adult5b) ## same as adult5
AIC(adult1b, adult1, adult1c, adult1d, adult1e, adult3, adult4, adult5, adult5b)
BIC(adult1b, adult1, adult1c, adult1d, adult1e, adult3, adult4, adult5, adult5b)
# adult 1b
resid <- resid(adult3)
fitted <- fitted(adult3)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)
## look pretty good

### looking at nonlinear
## generating 1st & 3rd Q for nonlinear models
summary(mort2.adult$Nadj)
summary(mort2.adult$Padj)
summary(mort2.adult$log10N)
summary(mort2.adult$log10P)

## simple random
mod7 <- mixmeta(lnRR ~ ns(Nadj, knots=c(0.5525,9.8925)) + 
                  ns(Padj, knots=c(0.500,2.500)), 
                data=mort2.adult,
                random= ~ 1| RefID, 
                method="ml", 
                control=list(addSlist=addS.adult))
summary(mod7) #
mod7a <- mixmeta(lnRR ~ ns(Nadj, knots=c(0.500)) + 
                  ns(Padj, knots=c(0.5485)), 
                data=mort2.adult,
                random= ~ 1| RefID, 
                method="ml", 
                control=list(addSlist=addS.adult))
summary(mod7a) #

mod7b <- mixmeta(lnRR ~ ns(log10N, knots=c(-0.2584 ,0.9953)) + 
                   ns(log10P, knots=c(-0.3010 ,0.3979 )), 
                 data=mort2.adult,
                 random= ~ 1| RefID, 
                 method="ml", 
                 control=list(addSlist=addS.adult))
summary(mod7b)
AIC(mod7a, mod7b)

mod7b2 <- mixmeta(lnRR ~ ns(log10N, knots=c(0.5485)) + 
                   ns(log10P, knots=c(-0.3010)), 
                 data=mort2.adult,
                 random= ~ 1| RefID, 
                 method="ml", 
                 control=list(addSlist=addS.adult))
summary(mod7b2)
AIC(mod7, mod7a, mod7b, mod7b2, adult3)
BIC(mod7, mod7a, mod7b, mod7b2, adult3)
summary(adult3)

### plotting with predictions
newdata.predictN.adult <- data.frame(log10N = mort2.adult$log10N,
                                     log10P = median(mort2.adult$log10P),
                                     Exposure.in.days = median(mort2.adult$Exposure.in.days))
newdata.predictP.adult <- data.frame(log10N = median(mort2.adult$log10N),
                                     log10P = mort2.adult$log10P,
                                     Exposure.in.days = median(mort2.adult$Exposure.in.days))

predN.adult <- predict(adult3, newdata=newdata.predictN.adult, ci=TRUE)
predP.adult <- predict(adult3, newdata=newdata.predictP.adult, ci=TRUE)
testCIN.adult <- cbind(mort2.adult, predN.adult)
testCIP.adult <- cbind (mort2.adult, predP.adult)
min_conc_N.adult <- testCIN.adult %>% 
  mutate(overlap0 = 0 >= ci.lb & 0 <= ci.ub) %>% 
  filter(overlap0==FALSE) %>% 
  summarize(min_N = min(Nadj),
            max_N = max(Nadj),
            min_P = min(Padj),
            max_P = max(Padj))

## acropora longicyathus was much lower than the others
## and now plotting for nitrogen as a predictor
ESNadult <- ggplot(testCIN.adult, 
                 aes(x = Nadj,
                     y = lnRR, color = Padj,
                     ymin = lnRR-lnSERR,
                     ymax = lnRR+lnSERR)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Log Risk Ratio (+/- s.e.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.15, 0.85),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) +  
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = Nadj, y = lnRR,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) 
ESNadult

ESPadult <- ggplot(testCIP.adult, 
                 aes(x = Padj,
                     y = lnRR, color = Nadj,
                     ymin = lnRR-lnSERR,
                     ymax = lnRR+lnSERR)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)")) +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = c(0.2, 0.85),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIN ("*mu*"M)")) +
  geom_ribbon(aes(x = Padj, y = lnRR,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Padj, y = fit), inherit.aes=FALSE) 
grid.arrange(ESNadult,ESPadult,ncol=2)
summary(adult3)

### quartile plots
adultplot3 <- mort2.adult%>%mutate(NadjQuantBins = cut(Nadj, 
                                                    breaks = unique(quantile(Nadj,probs=seq.int(0,1, by=1/4))), 
                                                    include.lowest=TRUE))
adultplot3 <- adultplot3%>%mutate(PadjQuantBins = cut(Padj, 
                                                  breaks = unique(quantile(Padj,probs=seq.int(0,1, by=1/4))), 
                                                  include.lowest=TRUE))
N <- ggplot(adultplot3, aes(x = NadjQuantBins, fill = lnRRpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Proportion of Reported Change in Adult Survival",
       fill = "") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(size = 15)) +
  scale_x_discrete(labels = c("[0.1,0.552]" = "0.1-0.6",
                              "(0.552,3.75]" = "0.6-3.8",
                              "(3.75,9.89]" = "3.8-9.9",
                              "(9.89,32.5]" = "9.9-32.5"))
N
P <- ggplot(adultplot3, aes(x = PadjQuantBins, fill = lnRRpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "",
       fill = "") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 15))  +
  scale_x_discrete(labels = c("[0.02,0.5]" = "0.02-0.5",
                              "(0.5,2.5]" = "0.5-2.5",
                              "(2.5,5]" = "2.5-5"))
P
grid.arrange(N,P,ncol=2)
