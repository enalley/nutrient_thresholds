#### NALLEY ET AL. 2022 NUTRIENT THRESHOLD ANALYSIS #####
### MAXIMUM PHOTOSYNTHETIC EFFICIENCY (MQY - Fv/Fm) ###

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
library(mgcv)


breaks_log10 <- function(x) {
  low <- floor(log10(min(x)))
  high <- ceiling(log10(max(x)))
  
  10^(seq.int(low, high))
}

getwd()
setwd("/")
allmqy <- as.data.frame(read.csv("ALL-MQY.csv", header=TRUE))
str(allmqy) ## everything is factor vs. numeric
allmqy$Coral.Type <- as.factor(allmqy$Coral.Type)
allmqy$Family<- as.factor(allmqy$Family)

# creating new variable that has both Genus and species
allmqy <- allmqy %>%
  mutate(Gsp = paste(Genus, Species, sep = " ")) 
head(allmqy)

### estimating ambient concentrations of N & P if listed as zero
allmqy$Nadj <- ifelse(allmqy$TOTAL.N == 0, 0.1, allmqy$TOTAL.N)
allmqy$Padj <- ifelse(allmqy$PO4 == 0, 0.02, allmqy$PO4)
allmqy$NO3adj <- ifelse(allmqy$NO3 == 0, 0.1, allmqy$NO3)
allmqy$NH4adj <- ifelse(allmqy$NH4 == 0, 0.1, allmqy$NH4)
### adding log variables
allmqy$log10N <- log10(allmqy$Nadj+1)
allmqy$log10P <- log10(allmqy$Padj+1)
allmqy$log10NO3 <- log10(allmqy$NO3adj+1)
allmqy$log10NH4 <- log10(allmqy$NH4adj+1)
### ordering data before calculating hedges
str(allmqy)
allmqy.ord <- allmqy[order(allmqy$Reference, allmqy$Control),]

### calculating hedges D 
covar <- by(allmqy.ord, allmqy.ord$Reference, function(x) 
  covar.smd(Response.level, SD, N, "smd", method="hedges", data = x))
allmqy.ord$smd <- unlist(lapply(covar, function(x) x$y))
allmqy.ord$vmd <- unlist(lapply(covar, function(x) x$v))
str(allmqy.ord) 
test <- allmqy.ord[,c("Reference", "Control", "smd")]

### removing controls for models
allmqy2 <- subset(allmqy.ord, Control=="exp")
summary(allmqy2)
allmqy2$smdpos <- as.factor(ifelse(allmqy2$smd > 0, "Increase", "Decrease"))

## 7 Studies with 12 References within
MQYTreatSummary <- allmqy2 %>% 
  summarize(Nmin = min(Nadj),
            Nmax = max(Nadj),
            Pmin = min(Padj),
            Pmax = max(Padj),
            DurMean = mean(Exposure.in.days),
            DurMin = min(Exposure.in.days),
            DurMax = max(Exposure.in.days))

### dotplot
ggplot(allmqy2, 
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
       size = "Treatment\nEffect Size", color = "Photosynthetic\nEfficiency") + 
  guides(color = guide_legend(order = 2, override.aes = list(size=4)),
         size = guide_legend(order = 1)) 

### COVARIANCE MATRICES
newlist <- list(NA)
for (i in seq(1,length(covar))) {
  newlist[i] <- list(covar[[i]]$S)
}

### LOOKING AT SMD
## N + P with simple random
mod1 <- mixmeta(smd ~ Nadj + Padj, 
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod1) 

mod1b <- mixmeta(smd ~ NO3adj + NH4adj + Padj, 
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod1b)

## log10N + log10P with simple random
mod2 <- mixmeta(smd ~ log10N + log10P, 
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod2) 

mod2b <- mixmeta(smd ~ log10NO3 + log10NH4 + log10P, 
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod2b) 

AIC(mod1, mod1b, mod2, mod2b)
BIC(mod1, mod1b, mod2, mod2b)
##  mod1 better

## Nadj with poly2 + Padj with poly2
mod3 <- mixmeta(smd ~ poly(Nadj, 2, raw = TRUE) + poly(Padj, 2, raw = TRUE), 
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod3) 

## Nadj with poly2 + Padj
mod4 <- mixmeta(smd ~ poly(Nadj, 2, raw = TRUE) + Padj, 
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod4) 

## Nadj with poly3 + Padj with poly3
mod5 <- mixmeta(smd ~ poly(Nadj, 3, raw = TRUE) + poly(Padj, 3, raw = TRUE),
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod5) 

## Nadj with poly3 + Padj 
mod6 <- mixmeta(smd ~ poly(Nadj, 3, raw = TRUE) + Padj, 
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod6) 

AIC(mod1, mod3, mod4, mod5, mod6) ## mod3 is best
BIC(mod1, mod3, mod4, mod5, mod6) ## mod3 is best

mod7 <- mixmeta(smd ~ poly(Nadj, 2, raw = TRUE) + poly(Padj, 2, raw = TRUE) + Exposure.in.days, 
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod7) 

allmqy2$log10Exp <- log10(allmqy2$Exposure.in.days)
mod9 <- mixmeta(smd ~ Nadj + Padj + Exposure.in.days, 
                 random = ~ 1 | Reference, 
                 data = allmqy2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod9) 


mod10 <- mixmeta(smd ~ Nadj + Padj + log10Exp, 
                 random = ~ 1 | Reference, 
                 data = allmqy2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod10) 

AIC(mod1, mod3, mod7, mod8, mod9, mod10) 
BIC(mod1, mod3, mod7, mod8, mod9, mod10)

summary(mod1)
### PLOTTING: PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS ####
newdata.predictN <- data.frame(Nadj = allmqy2$Nadj,
                               Padj = median(allmqy2$Padj))
newdata.predictP <- data.frame(Nadj = median(allmqy2$Nadj),
                               Padj = allmqy2$Padj)

predN <- predict(mod1, newdata=newdata.predictN, ci=TRUE)
predP <- predict(mod1, newdata=newdata.predictP, ci=TRUE)

testCIN <- cbind(allmqy2, predN)
testCIP <- cbind (allmqy2, predP)

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
        legend.position = c(0.15, 0.85),
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
        legend.position = c(0.15, 0.85),
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
EffP
grid.arrange(EffN,EffP,ncol=2)


#### colored by species for supp ####
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
        legend.position = c(0.2, 0.8),
        legend.text=element_text(size=7, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("")) 
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
EffP
grid.arrange(EffN,EffP,ncol=2)


#### colored by species for supp ####
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
        legend.position = c(0.2, 0.8),
        legend.text=element_text(size=9),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("")) 
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
  labs(color = expression("DIN ("*mu*"M)")) +
  geom_ribbon(aes(x = Padj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Padj, y = fit), inherit.aes=FALSE)
EffP
grid.arrange(EffN,EffP,ncol=2)

#### colored by species for supp ####
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
        legend.position = c(0.2, 0.8),
        legend.text=element_text(size=9),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("")) 
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
  labs(color = expression("DIN ("*mu*"M)")) +
  geom_ribbon(aes(x = Padj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Padj, y = fit), inherit.aes=FALSE)
EffP
grid.arrange(EffN,EffP,ncol=2)
