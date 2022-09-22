#### NALLEY ET AL. 2022 NUTRIENT THRESHOLD ANALYSIS #####
#### SYMB DENSITY (10^6 CELLS/CM^2) ####

rm(list= ls())
library(tidyverse)
library(mixmeta)
library(dosresmeta)
library(splines)
library(effsize)
library(lme4)
library(MuMIn)
library(RColorBrewer)
library(nlme)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(ggthemes)

breaks_log10 <- function(x) {
  low <- floor(log10(min(x)))
  high <- ceiling(log10(max(x)))
  
  10^(seq.int(low, high))
}

getwd()
setwd("/")
symb <- as.data.frame(read.csv("GitHub_ALL-CellsSA_7.9.22.csv", header=TRUE))
symb$Reference <- as.factor(symb$Reference)
symb$Study <- as.factor(symb$Study)
symb$Control <- as.factor(symb$Control)
symb$Genus <- as.factor(symb$Genus)
symb$Species <- as.factor(symb$Species)
symb$Ocean <- as.factor(symb$Ocean)
symb$pulse.press <- as.factor(symb$pulse.press)
symb$Coral.Type <- as.factor(symb$Coral.Type)
symb$Family<- as.factor(symb$Family)

str(symb) ## everything is correct as factor vs. numeric
symb <- symb[,c(1:17,20:21)]
symb <- as.data.frame(symb[!symb$Reference == "NU20",])

str(symb)
## only press experiments

# creating new variable that has both Genus and species
symb <- symb %>%
  mutate(Gsp = paste(Genus, Species, sep = " ")) 
head(symb)
### estimating ambient concentrations of N & P if listed as zero
symb$Nadj <- ifelse(symb$TOTAL.N == 0, 0.1, symb$TOTAL.N)
symb$Padj <- ifelse(symb$PO4 == 0, 0.02, symb$PO4)
symb$NO3adj <- ifelse(symb$NO3 == 0, 0.1, symb$NO3)
symb$NH4adj <- ifelse(symb$NH4 == 0, 0.1, symb$NH4)
### adding log variables
symb$log10N <- log10(symb$Nadj+1)
symb$log10P <- log10(symb$Padj+1)
symb$log10NO3 <- log10(symb$NO3adj+1)
symb$log10NH4 <- log10(symb$NH4adj+1)

### ordering data before calculating hedges
str(symb)
symb.ord <- symb[order(symb$Reference, symb$Control),]

### calculating hedges G (std. mean difference)
### Reference = experiment (within studies)
covar <- by(symb.ord, symb.ord$Reference, function(x) 
  covar.smd(Response.level, SD, N, "smd", method="hedges", data = x))
covar <- rlist::list.clean(covar)
symb.ord$smd <- unlist(lapply(covar, function(x) x$y))
symb.ord$vmd <- unlist(lapply(covar, function(x) x$v))
str(symb.ord) 
test <- symb.ord[,c("Reference", "Control", "smd")]

### removing controls for modeling
symb2 <- subset(symb.ord, Control=="exp")
str(symb2)
symb2$smdpos <- as.factor(ifelse(symb2$smd > 0, "Increase", "Decrease"))
str(symb2)
symb2$optdens <- as.factor(ifelse(symb2$Response.level > 3, "Exceeds", "Under"))
str(symb2)

## Summary of data 
str(symb2)
## 22 papers with 37 experiments within
SymbTreatSummary <- symb2 %>% 
  summarize(Nmin = min(Nadj),
            Nmax = max(Nadj),
            Pmin = min(Padj),
            Pmax = max(Padj),
            DurMean = mean(Exposure.in.days),
            DurMin = min(Exposure.in.days),
            DurMax = max(Exposure.in.days))

### dotplot
dotplotoverview <- ggplot(symb2, 
                          aes(x = Nadj, y = Padj, 
                              color = smdpos,
                              size = abs(smd),
                              shape = optdens)) + 
  ## surface water from station aloha 2019 annual report (surface ocean)
  annotate("point", x = 0.03, y = 0.03, colour = "black", size = 4, shape = 8) +
  ## ambient reef from Nakajima et al. 2015 (Malaysia)
  annotate("point", x = 0.75, y = 0.1, colour = "black", size = 4, shape = 8) +
  ## ambient reef from Silbiger et al. 2018 (Kaneohe Bay)
  annotate("point", x = .15, y = .15, color = "black", size = 4, shape = 8) + 
  ## ambient (but high) from AIMS reef in Fabricius et al. 2013 (GBR)
  annotate("point", x = .44, y = .38, color = "black", size = 4, shape = 8) + 
  ## high reef from Silbiger et al. 2018 (throughout Pacific)
  annotate("point", x = 7.6, y = 2.6, color = "black", size = 4, shape = 8) + 
  ## highest measured at SGD seep in Lubarsky et al. 2018 (Maunalua Bay)
  annotate("point", x = 32.39, y = 0.89, color = "black", size = 4, shape = 8) + 
  scale_x_log10(breaks = breaks_log10) +  
  scale_y_log10(breaks = breaks_log10) + annotation_logticks() +
  # geom_jitter(alpha = .5, position = position_jitter(width = .07, height = 0.03)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  scale_size(range = c(3,15)) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = expression("DIP ("*mu*"M)"),
       size = "Treatment\nEffect Size", color = "Symbiont Density", 
       shape = "Optimal Density") + 
  guides(shape = guide_legend(order = 3, override.aes = list(size=4)),
         color = guide_legend(order = 2, override.aes = list(size=4)),
         size = guide_legend(order = 1)) +
  scale_shape_manual(values=c(17, 16))
dotplotoverview

ggplot(symb2, 
       aes(x = Nadj, y = Padj, 
           color = smdpos,
           size = abs(smd),
           shape = optdens)) + 
  ## surface water from station aloha 2019 annual report (surface ocean)
  annotate("point", x = 0.03, y = 0.03, colour = "darkgray", size = 4, shape = 8) +
  ## ambient reef from Nakajima et al. 2015 (Malaysia)
  annotate("point", x = 0.75, y = 0.1, colour = "darkgray", size = 4, shape = 8) +
  ## ambient reef from Silbiger et al. 2018 (Kaneohe Bay)
  annotate("point", x = .15, y = .15, color = "darkgray", size = 4, shape = 8) + 
  ## ambient (but high) from AIMS reef in Fabricius et al. 2013 (GBR)
  annotate("point", x = .44, y = .38, color = "darkgray", size = 4, shape = 8) + 
  ## high reef from Silbiger et al. 2018 (throughout Pacific)
  annotate("point", x = 7.6, y = 2.6, color = "darkgray", size = 4, shape = 8) + 
  ## highest measured at SGD seep in Lubarsky et al. 2018 (Maunalua Bay)
  annotate("point", x = 32.39, y = 0.89, color = "darkgray", size = 4, shape = 8) + 
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
       size = "Treatment\nEffect Size", color = "Zooxanthellae\nDensity", 
       shape = "Optimal Density") + 
  guides(shape = guide_legend(order = 3, override.aes = list(size=4)),
         color = guide_legend(order = 2, override.aes = list(size=4)),
         size = guide_legend(order = 1)) +
  scale_shape_manual(values=c(17, 16))

### LOOKING AT SMD
### COVARIANCE MATRICES for meta-analysis
newlist <- list(NA)
for (i in seq(1,length(covar))) {
  newlist[i] <- list(covar[[i]]$S)
}

studyspcs <- symb2 %>% 
  group_by(Reference) %>% 
  count(Gsp)

## N + P with simple random
mod1 <- mixmeta(smd ~ Nadj + Padj, 
                random = ~ 1 | Reference, 
                data = symb2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod1) # 

## log10N + log10P with simple random
mod2 <- mixmeta(smd ~ log10N + log10P, 
                 random = ~ 1 | Reference, 
                 data = symb2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod2) # 

AIC(mod1, mod2)
BIC(mod1, mod2)


## log10N with poly2 + log10P with poly2
mod3 <- mixmeta(smd ~ poly(log10N, 2, raw = TRUE) + poly(log10P, 2, raw = TRUE), 
                      random = ~ 1 | Reference, 
                      data = symb2, method = "ml", 
                      control = list(addSlist = newlist))
summary(mod3) # 

## log10N with poly2 + log10P
mod4 <- mixmeta(smd ~ poly(log10N, 2, raw = TRUE) + log10P, 
                      random = ~ 1 | Reference, 
                      data = symb2, method = "ml", 
                      control = list(addSlist = newlist))
summary(mod4) # 

## log10N with poly3+ log10P with poly3
mod5 <- mixmeta(smd ~ poly(log10N, 3, raw = TRUE) + poly(log10P, 3, raw = TRUE),
                      random = ~ 1 | Reference, 
                      data = symb2, method = "ml", 
                      control = list(addSlist = newlist))
summary(mod5) #

## log10N with poly3 + log10P 
mod6 <- mixmeta(smd ~ poly(log10N, 3, raw = TRUE) + log10P, 
                     random = ~ 1 | Reference, 
                     data = symb2, method = "ml", 
                     control = list(addSlist = newlist))
summary(mod6) # 

AIC(mod2, mod3, mod4, mod5, mod6) 
BIC(mod2, mod3, mod4, mod5, mod6)

 
### log10N with poly2 + log10P + Gsp 
mod7 <- mixmeta(smd ~ poly(log10N, 3, raw = TRUE) + log10P + Gsp, 
                      random = ~ 1 | Reference, 
                      data = symb2, method = "ml", 
                      control = list(addSlist = newlist))
summary(mod7) 

## log10N with poly3 + log10P + Exposure in days
mod8 <- mixmeta(smd ~ poly(log10N, 3, raw = TRUE) + log10P + Exposure.in.days, 
                  random = ~ 1 | Reference, 
                  data = symb2, method = "ml", 
                  control = list(addSlist = newlist))
summary(mod8) 

AIC(mod6, mod7, mod8) 
BIC(mod6, mod7, mod8) 

##### looking at NO3 & NH4 separately #####
mod9 <- mixmeta(smd ~ log10NO3 + log10NH4 + log10P, 
                  random = ~ 1 | Reference, 
                  data = symb2, method = "ml", 
                  control = list(addSlist = newlist))
summary(mod9) 

mod10 <- mixmeta(smd ~ poly(log10NO3, 2, raw = TRUE) + 
                       poly(log10NH4, 2, raw = TRUE) + 
                       log10P, 
                     random = ~ 1 | Reference, 
                     data = symb2, method = "ml", 
                     control = list(addSlist = newlist))
summary(mod10) 

mod11 <- mixmeta(smd ~ poly(log10NO3, 3, raw = TRUE) + 
                   poly(log10NH4, 3, raw = TRUE) + 
                   log10P, 
                 random = ~ 1 | Reference, 
                 data = symb2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod11) 

mod12 <- mixmeta(smd ~ poly(log10NO3, 3, raw = TRUE) + 
                   poly(log10NH4, 2, raw = TRUE) + 
                   log10P, 
                 random = ~ 1 | Reference, 
                 data = symb2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod12) 

mod13 <- mixmeta(smd ~ poly(log10NO3, 3, raw = TRUE) + 
                   log10NH4 + 
                   log10P, 
                 random = ~ 1 | Reference, 
                 data = symb2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod13) 

mod14 <- mixmeta(smd ~ poly(log10NO3, 2, raw = TRUE) + 
                   log10NH4 + 
                   log10P, 
                 random = ~ 1 | Reference, 
                 data = symb2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod14) 

AIC(mod4, mod6, mod9, mod10, mod11, mod12, mod13, mod14) 
BIC(mod4, mod6, mod9, mod10, mod11, mod12, mod13, mod14) 

## and now testing with Gsp 
mod15 <- mixmeta(smd ~ poly(log10NO3, 2, raw = TRUE) + 
                   log10NH4 + 
                   log10P + Gsp, 
                  random = ~ 1 | Reference, 
                  data = symb2, method = "ml", 
                  control = list(addSlist = newlist))
summary(mod15) 

## and exposure duration
mod16 <- mixmeta(smd ~ poly(log10NO3, 2, raw = TRUE) + 
                   log10NH4 + 
                   log10P + Exposure.in.days, 
                  random = ~ 1 | Reference, 
                  data = symb2, method = "ml", 
                  control = list(addSlist = newlist))
summary(mod16) 
AIC(mod9, mod14, mod15, mod16)
BIC(mod9, mod14, mod15, mod16)

AIC(mod4, mod14)
BIC(mod4, mod14)

resid <- resid(mod14)
fitted <- fitted(mod14)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)

### PLOTTING: PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS ####
newdata.predictNO3 <- data.frame(log10NO3 = symb2$log10NO3,
                                 log10NH4 = median(symb2$log10NH4),
                                 log10P = median(symb2$log10P))
newdata.predictNH4 <- data.frame(log10NO3 = median(symb2$log10NO3),
                                 log10NH4 = symb2$log10NH4,
                                 log10P = median(symb2$log10P))
newdata.predictP <- data.frame(log10NO3 = median(symb2$log10NO3),
                               log10NH4 = median(symb2$log10NH4),
                               log10P = symb2$log10P)
predNO3 <- predict(mod14, newdata=newdata.predictNO3, ci=TRUE)
predNH4 <- predict(mod14, newdata=newdata.predictNH4, ci=TRUE)
predP <- predict(mod14, newdata=newdata.predictP, ci=TRUE)
testCINO3 <- cbind(symb2, predNO3)
testCINH4 <- cbind(symb2, predNH4)
testCIP <- cbind (symb2, predP)

## NO3 as a predictor
ESNO3lin <- ggplot(testCINO3, 
                   aes(x = NO3adj,
                       y = smd, color = Padj,
                       ymin = smd-vmd,
                       ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("Nitrate ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.85, 0.75),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) +  
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b")  +
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = NO3, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = NO3, y = fit), inherit.aes=FALSE)
ESNO3lin

## NH4 as a predictor
ESNH4lin <- ggplot(testCINH4, 
                   aes(x = NH4adj,
                       y = smd, color = Padj,
                       ymin = smd-vmd,
                       ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("Ammonium ("*mu*"M)"),
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
        legend.position = c(0.85, 0.75),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) +  
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b")  +
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = NH4, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = NH4, y = fit), inherit.aes=FALSE) 
ESNH4lin

## DIP as a predictor 
ESPlin <- ggplot(testCIP, 
                 aes(x = Padj,
                     y = smd, color = Nadj,
                     ymin = smd-vmd,
                     ymax = smd+vmd)) + 
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
        legend.position = c(0.15, 0.75),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIN ("*mu*"M)")) +
  geom_ribbon(aes(x = Padj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Padj, y = fit), inherit.aes=FALSE) 
ESPlin
grid.arrange(ESNO3lin,ESNH4lin, ESPlin,ncol=3)

#### and now colored by species for supplemental ####
## NO3 as a predictor
ESNO3lin2 <- ggplot(testCINO3, 
                   aes(x = NO3adj,
                       y = smd, color = Coral.Type,
                       ymin = smd-vmd,
                       ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("Nitrate ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.825, 0.725),
        legend.text=element_text(size=5, face = "italic"),
        legend.background = element_blank()) +  
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b")  +
  labs(color = element_blank()) +
  geom_ribbon(aes(x = NO3, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = NO3, y = fit), inherit.aes=FALSE)
ESNO3lin2

## NH4 as a predictor
ESNH4lin2 <- ggplot(testCINH4, 
                   aes(x = NH4adj,
                       y = smd, color = Coral.Type,
                       ymin = smd-vmd,
                       ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("Ammonium ("*mu*"M)"),
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
        legend.position = "none") +  
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b")  +
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = NH4, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = NH4, y = fit), inherit.aes=FALSE) 
ESNH4lin2

## DIP as a predictor 
ESPlin2 <- ggplot(testCIP, 
                 aes(x = Padj,
                     y = smd, color = Coral.Type,
                     ymin = smd-vmd,
                     ymax = smd+vmd)) + 
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
        legend.position = "none") +
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIN ("*mu*"M)")) +
  geom_ribbon(aes(x = Padj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Padj, y = fit), inherit.aes=FALSE) 
ESPlin2
grid.arrange(ESNO3lin2, ESNH4lin2, ESPlin2, ncol=3)

#### and now colored by Family for supplemental ####
## NO3 as a predictor
ESNO3lin2 <- ggplot(testCINO3, 
                    aes(x = NO3adj,
                        y = smd, color = Family,
                        ymin = smd-vmd,
                        ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("Nitrate ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.825, 0.825),
        legend.text=element_text(size=10, face = "italic"),
        legend.background = element_blank()) +  
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b")  +
  labs(color = element_blank()) +
  geom_ribbon(aes(x = NO3, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = NO3, y = fit), inherit.aes=FALSE)
ESNO3lin2

## NH4 as a predictor
ESNH4lin2 <- ggplot(testCINH4, 
                    aes(x = NH4adj,
                        y = smd, color = Family,
                        ymin = smd-vmd,
                        ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("Ammonium ("*mu*"M)"),
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
        legend.position = "none") +  
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b")  +
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = NH4, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = NH4, y = fit), inherit.aes=FALSE) 
ESNH4lin2

## DIP as a predictor 
ESPlin2 <- ggplot(testCIP, 
                  aes(x = Padj,
                      y = smd, color = Family,
                      ymin = smd-vmd,
                      ymax = smd+vmd)) + 
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
        legend.position = "none") +
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIN ("*mu*"M)")) +
  geom_ribbon(aes(x = Padj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Padj, y = fit), inherit.aes=FALSE) 
ESPlin2
grid.arrange(ESNO3lin2, ESNH4lin2, ESPlin2, ncol=3)

#### and now colored by branching/plating/mounding for supplemental ####
## NO3 as a predictor
ESNO3lin2 <- ggplot(testCINO3, 
                    aes(x = NO3adj,
                        y = smd, color = Coral.Type,
                        ymin = smd-vmd,
                        ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("Nitrate ("*mu*"M)"),
       y = "Effect size (Hedges' d +/- s.d.)") +
  geom_abline(intercept=0, slope=0) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.825, 0.85),
        legend.text=element_text(size=10, face = "italic"),
        legend.background = element_blank()) +  
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b")  +
  labs(color = element_blank()) +
  geom_ribbon(aes(x = NO3, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = NO3, y = fit), inherit.aes=FALSE)
ESNO3lin2

## NH4 as a predictor
ESNH4lin2 <- ggplot(testCINH4, 
                    aes(x = NH4adj,
                        y = smd, color = Coral.Type,
                        ymin = smd-vmd,
                        ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
  labs(x = expression("Ammonium ("*mu*"M)"),
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
        legend.position = "none") +  
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b")  +
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = NH4, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = NH4, y = fit), inherit.aes=FALSE) 
ESNH4lin2

## DIP as a predictor 
ESPlin2 <- ggplot(testCIP, 
                  aes(x = Padj,
                      y = smd, color = Coral.Type,
                      ymin = smd-vmd,
                      ymax = smd+vmd)) + 
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
        legend.position = "none") +
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIN ("*mu*"M)")) +
  geom_ribbon(aes(x = Padj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Padj, y = fit), inherit.aes=FALSE) 
ESPlin2
grid.arrange(ESNO3lin2, ESNH4lin2, ESPlin2, ncol=3)



