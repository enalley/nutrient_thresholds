#### NALLEY ET AL. 2022 NUTRIENT THRESHOLD ANALYSIS #####
#### ADULT SURVIVAL, LARVAL SURVIVAL/SETTLEMENT, & FERTILIZATION ####

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
getwd()
setwd("/")
mort <- as.data.frame(read.csv("GitHub_ALL-Mortality_updated3.28.22.csv", header=TRUE))
str(mort) ## everything is factor vs. numeric
mort$SD <- as.numeric(as.character(mort$SD))
mort$SE <- as.numeric(as.character(mort$SE))

# creating new variable that has both Genus and species
mort <- mort %>%
  mutate(Gsp = paste(Genus, Species, sep = " ")) 
head(mort)
mort$RRsig <- as.factor(ifelse(mort$RiskRatio > 1, 1, 0))
mort$lnRR <- log(mort$RiskRatio)

## will give errors if SERiskRatio = 0, so changing to 0.001
mort$SERiskRatio[mort$SERiskRatio == 0] <- 0.001 
mort$type <- "ir"
### adjusting for studies that don't list min value
mort$Nadj <- ifelse(mort$Total.N == 0, 0.1, mort$Total.N)
mort$Padj <- ifelse(mort$PO4 == 0, 0.02, mort$PO4)
### adding log variables
mort$log10N <- log10(mort$Nadj+1)
mort$log10P <- log10(mort$Padj+1)
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
  covar.logrr(cases=Cases, n=N, y=lnRR, v=I(SERiskRatio^2), type=type, data=x))
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


### generating regular effect sizes (SMD) for egg & larva
unique(mort.ord.larva$RefID)
covar.larva <- by(mort.ord.larva, mort.ord.larva$RefID, function(x) 
  covar.smd(Response, SD, N, measure = "smd", method="hedges", data = x))
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


##### adult ####
adult1 <- mixmeta(lnRR ~  log10N + log10P,
                  random = ~ 1 | RefID,
                  data=mort2.adult, method="ml",
                  control=list(addSlist=addS.adult))
summary(adult1) 

adult2 <- mixmeta(lnRR ~ Nadj + Padj,
                   random = ~ 1 | RefID,
                   data=mort2.adult, method="ml",
                   control=list(addSlist=addS.adult))
summary(adult2) 
AIC(adult1, adult2)
BIC(adult1, adult2)

## log10N with poly2 + log10P with poly2
adult3 <- mixmeta(lnRR ~ poly(log10N, 2, raw = TRUE) + 
                    poly(log10P, 2, raw = TRUE), 
                random = ~ 1 | RefID, 
                data = mort2.adult, method = "ml", 
                control = list(addSlist = addS.adult))
summary(adult3)

## log10N with poly2 + log10P
adult4 <- mixmeta(lnRR ~ poly(log10N, 2, raw = TRUE) + log10P, 
                random = ~ 1 | RefID, 
                data = mort2.adult, method = "ml", 
                control = list(addSlist = addS.adult))
summary(adult4)

## log10N with poly3+ log10P with poly3
adult5 <- mixmeta(lnRR ~ poly(log10N, 3, raw = TRUE) + 
                    poly(log10P, 3, raw = TRUE),
                random = ~ 1 | RefID, 
                data = mort2.adult, method = "ml", 
                control = list(addSlist = addS.adult))
summary(adult5)

## log10N with poly3 + log10P 
adult6 <- mixmeta(lnRR ~ poly(log10N, 3, raw = TRUE) + log10P, 
                random = ~ 1 | RefID, 
                data = mort2.adult, method = "ml", 
                control = list(addSlist = addS.adult))
summary(adult6) 

AIC(adult1, adult3, adult4, adult5, adult6) 
BIC(adult1, adult3, adult4, adult5, adult6)

adult7 <- mixmeta(lnRR ~ log10N + log10P + Gsp, 
                random = ~ 1 | RefID, 
                data = mort2.adult, method = "ml", 
                control = list(addSlist = addS.adult))
summary(adult7) 

mort2.adult$log10Exp <- log10(mort2.adult$Exposure.in.days+1)
adult8 <- mixmeta(lnRR ~ log10N + log10P + Exposure.in.days, 
                 random = ~ 1 | RefID, 
                 data = mort2.adult, method = "ml", 
                 control = list(addSlist = addS.adult))
summary(adult8) 

adult9 <- mixmeta(lnRR ~ log10N + log10P + log10Exp, 
                 random = ~ 1 | RefID, 
                 data = mort2.adult, method = "ml", 
                 control = list(addSlist = addS.adult))
summary(adult9) 

adult10 <- mixmeta(lnRR ~ log10N + log10P + Exposure.in.days +
                     log10N*log10P + log10P*Exposure.in.days +
                     log10N*Exposure.in.days, 
                  random = ~ 1 | RefID, 
                  data = mort2.adult, method = "ml", 
                  control = list(addSlist = addS.adult))
summary(adult10)

adult11 <- mixmeta(lnRR ~ log10N*Exposure.in.days + 
                     log10P*Exposure.in.days, 
                   random = ~ 1 | RefID, 
                   data = mort2.adult, method = "ml", 
                   control = list(addSlist = addS.adult))
summary(adult11)

adult12 <- mixmeta(lnRR ~ log10N*log10Exp + 
                     log10P*log10Exp, 
                   random = ~ 1 | RefID, 
                   data = mort2.adult, method = "ml", 
                   control = list(addSlist = addS.adult))
summary(adult12)

AIC(adult1, adult4, adult7, adult8, adult9, adult10, adult11, adult12)
BIC(adult1, adult4, adult7, adult8, adult9, adult10, adult11, adult12)

### PLOTTING: PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS ####
newdata.predictN <- data.frame(log10N = mort2.adult$log10N,
                               log10P = median(mort2.adult$log10P),
                               Exposure.in.days = median(mort2.adult$Exposure.in.days))
newdata.predictP <- data.frame(log10N = median(mort2.adult$log10N),
                               log10P = mort2.adult$log10P,
                               Exposure.in.days = median(mort2.adult$Exposure.in.days))
newdata.predictE <- data.frame(log10N = median(mort2.adult$log10N),
                               log10P = median(mort2.adult$log10P),
                               Exposure.in.days = mort2.adult$Exposure.in.days)
predN <- predict(adult8, newdata=newdata.predictN, ci=TRUE)
predP <- predict(adult8, newdata=newdata.predictP, ci=TRUE)
predE <- predict(adult8, newdata=newdata.predictE, ci=TRUE)
testCIN <- cbind(mort2.adult, predN)
testCIP <- cbind (mort2.adult, predP)
testCIE <- cbind (mort2.adult, predE)

EffN <- ggplot(testCIN, 
               aes(x = Nadj,
                   y = lnRR, color = Padj,
                   ymin = lnRR-SERiskRatio,
                   ymax = lnRR+SERiskRatio)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Log Risk Ratio (+/- s.e.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.2, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIP ("*mu*"M)")) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP, 
               aes(x = Padj,
                   y = lnRR, color = Nadj,
                   ymin = lnRR-SERiskRatio,
                   ymax = lnRR+SERiskRatio)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Log Risk Ratio (+/- s.e.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = c(0.2, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIN ("*mu*"M)")) 
EffP

## and Exp as a predictor
EffE <- ggplot(testCIE, 
               aes(x = Exposure.in.days,
                   y = lnRR, color = Padj,
                   ymin = lnRR-SERiskRatio,
                   ymax = lnRR+SERiskRatio)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("Exposure (days)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = c(0.2, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = Exposure.in.days, y = lnRR,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Exposure.in.days, y = fit), inherit.aes=FALSE)
EffE
grid.arrange(EffN,EffP,EffE,ncol=3)


### COLORED BY SPECIES FOR SUPPLEMENTAL
EffN <- ggplot(testCIN, 
               aes(x = Nadj,
                   y = lnRR, color = Gsp,
                   ymin = lnRR-SERiskRatio,
                   ymax = lnRR+SERiskRatio)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Log Risk Ratio (+/- s.e.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        legend.text=element_text(size=12, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("Species")) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP, 
               aes(x = Padj,
                   y = lnRR, color = Gsp,
                   ymin = lnRR-SERiskRatio,
                   ymax = lnRR+SERiskRatio)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.text=element_text(size=12, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression(""))  
EffP

EffE <- ggplot(testCIE, 
               aes(x = Exposure.in.days,
                   y = lnRR, color = Gsp,
                   ymin = lnRR-SERiskRatio,
                   ymax = lnRR+SERiskRatio)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("Exposure (days)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = c(0.3, 0.3),
        legend.text=element_text(size=12, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  labs(color = expression("")) +
  geom_ribbon(aes(x = Exposure.in.days, y = lnRR,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Exposure.in.days, y = fit), inherit.aes=FALSE)
EffE
grid.arrange(EffN,EffP,EffE,ncol=3)

### COLORED BY Family FOR SUPPLEMENTAL
EffN <- ggplot(testCIN, 
               aes(x = Nadj,
                   y = lnRR, color = Family,
                   ymin = lnRR-SERiskRatio,
                   ymax = lnRR+SERiskRatio)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Log Risk Ratio (+/- s.e.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        legend.text=element_text(size=12, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("Species")) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP, 
               aes(x = Padj,
                   y = lnRR, color = Family,
                   ymin = lnRR-SERiskRatio,
                   ymax = lnRR+SERiskRatio)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.text=element_text(size=12, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression(""))  
EffP

EffE <- ggplot(testCIE, 
               aes(x = Exposure.in.days,
                   y = lnRR, color = Family,
                   ymin = lnRR-SERiskRatio,
                   ymax = lnRR+SERiskRatio)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("Exposure (days)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = c(0.3, 0.3),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  labs(color = expression("")) +
  geom_ribbon(aes(x = Exposure.in.days, y = lnRR,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Exposure.in.days, y = fit), inherit.aes=FALSE)
EffE
grid.arrange(EffN,EffP,EffE,ncol=3)

### COLORED BY Morphology FOR SUPPLEMENTAL
EffN <- ggplot(testCIN,  
               aes(x = Nadj,
                   y = lnRR, color = Coral.Type,
                   ymin = lnRR-SERiskRatio,
                   ymax = lnRR+SERiskRatio)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Log Risk Ratio (+/- s.e.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        legend.text=element_text(size=12, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("Species")) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP, 
               aes(x = Padj,
                   y = lnRR, color = Coral.Type,
                   ymin = lnRR-SERiskRatio,
                   ymax = lnRR+SERiskRatio)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.text=element_text(size=12, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression(""))  
EffP

EffE <- ggplot(testCIE, 
               aes(x = Exposure.in.days,
                   y = lnRR, color = Coral.Type,
                   ymin = lnRR-SERiskRatio,
                   ymax = lnRR+SERiskRatio)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("Exposure (days)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = c(0.3, 0.3),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  labs(color = expression("")) +
  geom_ribbon(aes(x = Exposure.in.days, y = lnRR,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Exposure.in.days, y = fit), inherit.aes=FALSE)
EffE
grid.arrange(EffN,EffP,EffE,ncol=3)


###### LARVA ####
## choosing to use smd vs. lnRR for consistency with other analyses
mort2b.larva.test <- mort2b.larva[!(mort2b.larva$RefID=="NU11d" | mort2b.larva$RefID=="NU73e"),]
mort2b.larva.test
newlist.larva.test <- newlist.larva[-c(5,17)]

## N + P with simple random
larva1 <- mixmeta(smd ~ Nadj + Padj, 
                  random = ~ 1 | RefID, 
                  data = mort2b.larva, method = "ml", 
                  control = list(addSlist = newlist.larva))
summary(larva1) # 

larva1b <- mixmeta(smd ~ Nadj + Padj, 
                  random = ~ 1 | RefID, 
                  data = mort2b.larva.test, method = "ml", 
                  control = list(addSlist = newlist.larva.test))
summary(larva1b) # 

larva1c <- mixmeta(lnRR ~ Nadj + Padj, 
                  random = ~ 1 | RefID, 
                  data = mort2b.larva, method = "ml", 
                  control = list(addSlist = addS.larva))
summary(larva1c)

## log10N + log10P with simple random
larva2 <- mixmeta(smd ~ log10N + log10P, 
                  random = ~ 1 | RefID, 
                  data = mort2b.larva, method = "ml", 
                  control = list(addSlist = newlist.larva))
summary(larva2)

larva2b <- mixmeta(smd ~ log10N + log10P, 
                random = ~ 1 | RefID, 
                data = mort2b.larva.test, method = "ml", 
                control = list(addSlist = newlist.larva.test))
summary(larva2b) ## 


AIC(larva1, larva2)
BIC(larva1, larva2)

## log10N with poly2 + log10P with poly2
larva3 <- mixmeta(smd ~ poly(log10N, 2, raw = TRUE) + 
                    poly(log10P, 2, raw = TRUE), 
                  random = ~ 1 | RefID, 
                  data = mort2b.larva, method = "ml", 
                  control = list(addSlist = newlist.larva))
summary(larva3) # 

## log10N with poly2 + log10P
larva4 <- mixmeta(smd ~ poly(log10N, 2, raw = TRUE) + log10P, 
                  random = ~ 1 | RefID, 
                  data = mort2b.larva, method = "ml", 
                  control = list(addSlist = newlist.larva))
summary(larva4) # 

## log10N with poly3+ log10P with poly3
larva5 <- mixmeta(smd ~ poly(log10N, 3, raw = TRUE) + 
                    poly(log10P, 3, raw = TRUE),
                  random = ~ 1 | RefID, 
                  data = mort2b.larva, method = "ml", 
                  control = list(addSlist = newlist.larva))
summary(larva5) # 

## log10N with poly3 + log10P 
larva6 <- mixmeta(smd ~ poly(log10N, 3, raw = TRUE) + log10P, 
                  random = ~ 1 | RefID, 
                  data = mort2b.larva, method = "ml", 
                  control = list(addSlist = newlist.larva))
summary(larva6) # 

AIC(larva1, larva2, larva3, larva4, larva5, larva6) 
BIC(larva1, larva2, larva3, larva4, larva5, larva6)

larva7 <- mixmeta(smd ~ log10N + log10P + Gsp, 
                  random = ~ 1 | RefID, 
                  data = mort2b.larva, method = "ml", 
                  control = list(addSlist = newlist.larva))
summary(larva7) ## 

mort2b.larva$log10Exp <- log10(mort2b.larva$Exposure.in.days+1)
larva8 <- mixmeta(smd ~ log10N + log10P + Exposure.in.days, 
                  random = ~ 1 | RefID, 
                  data = mort2b.larva, method = "ml", 
                  control = list(addSlist = newlist.larva))
summary(larva8) ## 

larva9 <- mixmeta(smd ~ log10N + log10P + log10Exp, 
                  random = ~ 1 | RefID, 
                  data = mort2b.larva, method = "ml", 
                  control = list(addSlist = newlist.larva))
summary(adult9)

AIC(larva1, larva7, larva8, larva9)
BIC(larva1, larva7, larva8, larva9)

### PLOTTING: PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS ####
newdata.predictN <- data.frame(Nadj = mort2b.larva$Nadj,
                               Padj = median(mort2b.larva$Padj))
newdata.predictP <- data.frame(Nadj = median(mort2b.larva$Nadj),
                               Padj = mort2b.larva$Padj)
predN <- predict(larva1, newdata=newdata.predictN, ci=TRUE)
predP <- predict(larva1, newdata=newdata.predictP, ci=TRUE)
testCIN <- cbind(mort2b.larva, predN)
testCIP <- cbind (mort2b.larva, predP)

EffN <- ggplot(testCIN,
               aes(x = Nadj,
                   y = smd, color = Padj,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.15, 0.8),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = Nadj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP,
               aes(x = Padj,
                   y = smd, color = Nadj,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position =  c(0.85, 0.8),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIN ("*mu*"M)")) 
EffP
grid.arrange(EffN,EffP,ncol=2)

### COLORED BY SPECIES FOR SUPPLEMENTAL
EffN <- ggplot(testCIN[!(testCIN$RefID=="NU11d" | testCIN$RefID=="NU73e"),], 
               aes(x = Nadj,
                   y = smd, color = Gsp,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        legend.text=element_text(size=12, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("Species")) +
  geom_ribbon(aes(x = Nadj, y = lnRR,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP[!(testCIN$RefID=="NU11d" | testCIN$RefID=="NU73e"),], 
               aes(x = Padj,
                   y = smd, color = Gsp,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = c(0.75, 0.8),
        legend.text=element_text(size=12, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("Species"))  
EffP
grid.arrange(EffN,EffP,ncol=2)


### COLORED BY tax family FOR SUPPLEMENTAL
EffN <- ggplot(testCIN[!(testCIN$RefID=="NU11d" | testCIN$RefID=="NU73e"),], 
               aes(x = Nadj,
                   y = smd, color = Family,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        legend.text=element_text(size=12, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("Species")) +
  geom_ribbon(aes(x = Nadj, y = lnRR,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP[!(testCIN$RefID=="NU11d" | testCIN$RefID=="NU73e"),], 
               aes(x = Padj,
                   y = smd, color = Family,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = c(0.75, 0.8),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression(""))  
EffP
grid.arrange(EffN,EffP,ncol=2)

### COLORED BY morph FOR SUPPLEMENTAL
EffN <- ggplot(testCIN[!(testCIN$RefID=="NU11d" | testCIN$RefID=="NU73e"),], 
               aes(x = Nadj,
                   y = smd, color = Coral.Type,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        legend.text=element_text(size=12, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("Species")) +
  geom_ribbon(aes(x = Nadj, y = lnRR,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP[!(testCIN$RefID=="NU11d" | testCIN$RefID=="NU73e"),], 
               aes(x = Padj,
                   y = smd, color = Coral.Type,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = c(0.75, 0.8),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression(""))  
EffP
grid.arrange(EffN,EffP,ncol=2)

###### EGG ####
## N + P with simple random
egg1 <- mixmeta(smd ~ Nadj + Padj, 
                  random = ~ 1 | RefID, 
                  data = mort2b.egg, method = "ml", 
                  control = list(addSlist = newlist.egg))
summary(egg1) # 

egg1b <- mixmeta(lnRR ~ Nadj + Padj, 
                   random = ~ 1 | RefID, 
                   data = mort2b.egg, method = "ml", 
                   control = list(addSlist = addS.egg))
summary(egg1b)

## log10N + log10P with simple random
egg2 <- mixmeta(smd ~ log10N + log10P, 
                  random = ~ 1 | RefID, 
                  data = mort2b.egg, method = "ml", 
                  control = list(addSlist = newlist.egg))
summary(egg2)

egg2b <- mixmeta(lnRR ~ log10N + log10P, 
                 random = ~ 1 | RefID, 
                 data = mort2b.egg, method = "ml", 
                 control = list(addSlist = addS.egg))
summary(egg2b)

AIC(egg1, egg2)
BIC(egg1, egg2)

## log10N with poly2 + log10P with poly2
egg3 <- mixmeta(smd ~ poly(Nadj, 2, raw = TRUE) + 
                    poly(Padj, 2, raw = TRUE), 
                  random = ~ 1 | RefID, 
                  data = mort2b.egg, method = "ml", 
                  control = list(addSlist = newlist.egg))
summary(egg3) # 

## log10N with poly2 + log10P
egg4 <- mixmeta(smd ~ poly(Nadj, 2, raw = TRUE) + Padj, 
                  random = ~ 1 | RefID, 
                  data = mort2b.egg, method = "ml", 
                  control = list(addSlist = newlist.egg))
summary(egg4) # 

egg4b <- mixmeta(smd ~ poly(Padj, 2, raw = TRUE) + Nadj, 
                random = ~ 1 | RefID, 
                data = mort2b.egg, method = "ml", 
                control = list(addSlist = newlist.egg))
summary(egg4b) 

## log10N with poly3+ log10P with poly3
egg5 <- mixmeta(smd ~ poly(Nadj, 3, raw = TRUE) + 
                    poly(Padj, 3, raw = TRUE),
                  random = ~ 1 | RefID, 
                  data = mort2b.egg, method = "ml", 
                  control = list(addSlist = newlist.egg))
summary(egg5) # 

## log10N with poly3 + log10P 
egg6 <- mixmeta(smd ~ poly(Nadj, 3, raw = TRUE) + Padj, 
                  random = ~ 1 | RefID, 
                  data = mort2b.egg, method = "ml", 
                  control = list(addSlist = newlist.egg))
summary(egg6) # 

egg6b <- mixmeta(smd ~ poly(Padj, 3, raw = TRUE) + Nadj, 
                random = ~ 1 | RefID, 
                data = mort2b.egg, method = "ml", 
                control = list(addSlist = newlist.egg))
summary(egg6b)

AIC(egg1, egg2, egg3, egg4, egg4b, egg5, egg6, egg6b) 
BIC(egg1, egg2, egg3, egg4, egg4b, egg5, egg6, egg6b)

egg7 <- mixmeta(smd ~ log10N + log10P + Gsp, 
                  random = ~ 1 | RefID, 
                  data = mort2b.egg, method = "ml", 
                  control = list(addSlist = newlist.egg))
summary(egg7) ## 

egg8 <- mixmeta(smd ~ log10N + log10P + Exposure.in.days, 
                  random = ~ 1 | RefID, 
                  data = mort2b.egg, method = "ml", 
                  control = list(addSlist = newlist.egg))
summary(egg8) ## 

AIC(egg1, egg7, egg8)
BIC(egg1, egg7, egg8)

resid <- resid(egg1)
fitted <- fitted(egg1)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)

### PLOTTING: PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS ####
newdata.predictN <- data.frame(Nadj = mort2b.egg$Nadj,
                               Padj = median(mort2b.egg$Padj))
newdata.predictP <- data.frame(Nadj = median(mort2b.egg$Nadj),
                               Padj = mort2b.egg$log10P)
predN <- predict(egg1, newdata=newdata.predictN, ci=TRUE)
predP <- predict(egg1, newdata=newdata.predictP, ci=TRUE)
testCIN <- cbind(mort2b.egg, predN)
testCIP <- cbind (mort2b.egg, predP)

EffN <- ggplot(testCIN,
               aes(x = Nadj,
                   y = smd, color = Padj,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.75, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = Nadj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP,
               aes(x = Padj,
                   y = smd, color = Nadj,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position =  c(0.15, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIN ("*mu*"M)")) 
EffP
grid.arrange(EffN,EffP,ncol=2)

### COLORED BY SPECIES FOR SUPPLEMENTAL
EffN <- ggplot(testCIN, 
               aes(x = Nadj,
                   y = smd, color = Gsp,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.75, 0.2),
        legend.text=element_text(size=12, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("Species")) +
  geom_ribbon(aes(x = Nadj, y = lnRR,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) 

EffN

## and P as a predictor
EffP <- ggplot(testCIP, 
               aes(x = Padj,
                   y = smd, color = Gsp,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.text=element_text(size=12, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("Species"))  
EffP
grid.arrange(EffN,EffP,ncol=2)


### COLORED BY tax family FOR SUPPLEMENTAL
EffN <- ggplot(testCIN, 
               aes(x = Nadj,
                   y = smd, color = Family,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.75, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("")) +
  geom_ribbon(aes(x = Nadj, y = lnRR,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP, 
               aes(x = Padj,
                   y = smd, color = Family,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression(""))  
EffP
grid.arrange(EffN,EffP,ncol=2)

### COLORED BY morph FOR SUPPLEMENTAL
EffN <- ggplot(testCIN, 
               aes(x = Nadj,
                   y = smd, color = Coral.Type,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.75, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("")) +
  geom_ribbon(aes(x = Nadj, y = lnRR,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP, 
               aes(x = Padj,
                   y = smd, color = Coral.Type,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression(""))  
EffP
grid.arrange(EffN,EffP,ncol=2)

