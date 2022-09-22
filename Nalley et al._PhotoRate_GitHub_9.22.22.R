#### NALLEY ET AL. 2022 NUTRIENT THRESHOLD ANALYSIS #####
### Photosynthetic Rate ###

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

breaks_log10 <- function(x) {
  low <- floor(log10(min(x)))
  high <- ceiling(log10(max(x)))
  
  10^(seq.int(low, high))
}

getwd()
setwd("/")
allpho <- as.data.frame(read.csv("ALL-PhotoSA.csv", header=TRUE))
str(allpho) ## everything is factor vs. numeric
allpho$Coral.Type <- as.factor(allpho$Coral.Type)
allpho$Family<- as.factor(allpho$Family)
allpho <- allpho[,-c(21:27)]

# creating new variable that has both Genus and species
allpho <- allpho %>%
  mutate(Gsp = paste(Genus, Species, sep = " ")) 
head(allpho)

### estimating ambient concentrations of N & P if listed as zero
allpho$Nadj <- ifelse(allpho$TOTAL.N == 0, 0.1, allpho$TOTAL.N)
allpho$Padj <- ifelse(allpho$PO4 == 0, 0.02, allpho$PO4)
allpho$NO3adj <- ifelse(allpho$NO3 == 0, 0.1, allpho$NO3)
allpho$NH4adj <- ifelse(allpho$NH4 == 0, 0.1, allpho$NH4)
### adding log variables
allpho$log10N <- log10(allpho$Nadj+1)
allpho$log10P <- log10(allpho$Padj+1)
allpho$log10NO3 <- log10(allpho$NO3adj+1)
allpho$log10NH4 <- log10(allpho$NH4adj+1)
### ordering data before calculating hedges
str(allpho)
allpho.ord <- allpho[order(allpho$Reference, allpho$Control),]

### calculating hedges D 
covar <- by(allpho.ord, allpho.ord$Reference, function(x) 
  covar.smd(Response.level, SD, N, "smd", method="hedges", data = x))
allpho.ord$smd <- unlist(lapply(covar, function(x) x$y))
allpho.ord$vmd <- unlist(lapply(covar, function(x) x$v))
str(allpho.ord) 
test <- allpho.ord[,c("Reference", "Control", "smd")]

### removing controls for models
allpho2 <- subset(allpho.ord, Control=="exp")
summary(allpho2)
allpho2$smdpos <- as.factor(ifelse(allpho2$smd > 0, "Increase", "Decrease"))
str(allpho2)

## 11 Studies with 9 References within
PhotoTreatSummary <- allpho2 %>% 
  summarize(Nmin = min(Nadj),
            Nmax = max(Nadj),
            Pmin = min(Padj),
            Pmax = max(Padj),
            DurMean = mean(Exposure.in.days),
            DurMin = min(Exposure.in.days),
            DurMax = max(Exposure.in.days))

### dotplot
ggplot(allpho2, 
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
       size = "Treatment\nEffect Size", color = "Photo. Rate") +
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
                data = allpho2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod1) 

## log10N + log10P with simple random
mod2 <- mixmeta(smd ~ log10N + log10P, 
                random = ~ 1 | Reference, 
                data = allpho2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod2) 

AIC(mod1, mod2)
BIC(mod1, mod2)

## log10N with poly2 + log10P with poly2
mod3 <- mixmeta(smd ~ poly(log10N, 2, raw = TRUE) + poly(log10P, 2, raw = TRUE), 
                random = ~ 1 | Reference, 
                data = allpho2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod3)

## log10N with poly2 + log10P
mod4 <- mixmeta(smd ~ poly(log10N, 2, raw = TRUE) + log10P, 
                random = ~ 1 | Reference, 
                data = allpho2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod4) 

## log10N with poly3+ log10P with poly3
mod5 <- mixmeta(smd ~ poly(log10N, 3, raw = TRUE) + poly(log10P, 3, raw = TRUE),
                random = ~ 1 | Reference, 
                data = allpho2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod5) 

## log10N with poly3 + log10P 
mod6 <- mixmeta(smd ~ poly(log10N, 3, raw = TRUE) + log10P, 
                random = ~ 1 | Reference, 
                data = allpho2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod6)

AIC(mod2, mod3, mod4, mod5, mod6) 
BIC(mod2, mod3, mod4, mod5, mod6) 

### log10N with poly2 + log10P + Gsp 
mod7 <- mixmeta(smd ~ poly(log10N, 2, raw = TRUE) + log10P + Gsp, 
                random = ~ 1 | Reference, 
                data = allpho2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod7)  

mod8 <- mixmeta(smd ~ poly(log10N, 2, raw = TRUE) + log10P + Exposure.in.days, 
                random = ~ 1 | Reference, 
                data = allpho2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod8) 

## log10N with poly3 + log10P + Exposure in days
mod9 <- mixmeta(smd ~ log10N + log10P + Gsp, 
                random = ~ 1 | Reference, 
                data = allpho2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod9) ## 

allpho2$log10Exp <- log10(allpho2$Exposure.in.days+1)
mod10 <- mixmeta(smd ~ log10N + log10P + Exposure.in.days, 
                 random = ~ 1 | Reference, 
                 data = allpho2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod10) ## 


mod11 <- mixmeta(smd ~ log10N + log10P + log10Exp, 
                 random = ~ 1 | Reference, 
                 data = allpho2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod11) # 

AIC(mod2, mod4, mod7, mod8, mod9, mod10, mod11) 
BIC(mod2, mod4, mod7, mod8, mod9, mod10, mod11)
summary(mod7)
summary(mod4)

## with NO3 & NH4 sep
mod12 <- mixmeta(smd ~ log10NO3 + log10NH4 + log10P, 
                 random = ~ 1 | Reference, 
                 data = allpho2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod12) 

mod13 <- mixmeta(smd ~ poly(log10NO3, 2, raw = TRUE) + 
                   poly(log10NH4, 2, raw = TRUE) + 
                   log10P, 
                 random = ~ 1 | Reference, 
                 data = allpho2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod13) 

mod14 <- mixmeta(smd ~ poly(log10NO3, 2, raw = TRUE) + 
                   log10NH4 + 
                   log10P, 
                 random = ~ 1 | Reference, 
                 data = allpho2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod14) 

mod15 <- mixmeta(smd ~ log10NO3 + log10NH4 + log10P + Gsp, 
                 random = ~ 1 | Reference, 
                 data = allpho2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod15) 

mod16 <- mixmeta(smd ~ log10NO3 + log10NH4 + log10P + Exposure.in.days, 
                 random = ~ 1 | Reference, 
                 data = allpho2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod16) 

AIC(mod4, mod7, mod12, mod13, mod14, mod15, mod16) 
BIC(mod4, mod7, mod12, mod13, mod14, mod15, mod16) 

### PLOTTING: PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS ####
newdata.predictNO3 <- data.frame(log10NO3 = allpho2$log10NO3,
                               log10NH4 = median(allpho2$log10NH4),
                               log10P = median(allpho2$log10P))
newdata.predictNH4 <- data.frame(log10NO3 = median(allpho2$log10NO3),
                                 log10NH4 = allpho2$log10NH4,
                                 log10P = median(allpho2$log10P))
newdata.predictP <- data.frame(log10NO3 = median(allpho2$log10NO3),
                               log10NH4 = median(allpho2$log10NH4),
                               log10P = allpho2$log10P)
predNO3 <- predict(mod12, newdata=newdata.predictNO3, ci=TRUE)
predNH4 <- predict(mod12, newdata=newdata.predictNH4, ci=TRUE)
predP <- predict(mod12, newdata=newdata.predictP, ci=TRUE)
testCINO3 <- cbind(allpho2, predNO3)
testCINH4 <- cbind(allpho2, predNH4)
testCIP <- cbind (allpho2, predP)

EffNO3 <- ggplot(testCINO3, 
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
        legend.position = c(0.2, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = NO3adj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = NO3adj, y = fit), inherit.aes=FALSE) 
EffNO3

EffNH4 <- ggplot(testCINH4, 
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
        legend.position = c(0.15, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIP ("*mu*"M)")) 
EffNH4

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
        legend.position = c(0.85, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIN ("*mu*"M)")) 
EffP
grid.arrange(EffNO3,EffNH4,EffP,ncol=3)
summary(mod12)


#### and colored by species for supplemental ####
EffNO3 <- ggplot(testCINO3, 
                 aes(x = NO3adj,
                     y = smd, color = Gsp,
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
        legend.position = c(0.75, 0.2),
        legend.text=element_text(size=10, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("")) +
  geom_ribbon(aes(x = NO3adj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = NO3adj, y = fit), inherit.aes=FALSE) 
EffNO3

EffNH4 <- ggplot(testCINH4, 
                 aes(x = NH4adj,
                     y = smd, color = Gsp,
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
        legend.position = "none",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIP ("*mu*"M)")) 
EffNH4

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
  labs(color = expression("Species")) 
EffP
grid.arrange(EffNO3,EffNH4,EffP,ncol=3)

#### col by tax family for supplemental ####
EffNO3 <- ggplot(testCINO3, 
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
        legend.position = c(0.75, 0.2),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("")) +
  geom_ribbon(aes(x = NO3adj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = NO3adj, y = fit), inherit.aes=FALSE) 
EffNO3

EffNH4 <- ggplot(testCINH4, 
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
        legend.position = "none",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIP ("*mu*"M)")) 
EffNH4

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
  labs(color = expression("Species")) 
EffP
grid.arrange(EffNO3,EffNH4,EffP,ncol=3)


#### col by coral type for supplemental ####
EffNO3 <- ggplot(testCINO3, 
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
        legend.position = c(0.75, 0.2),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("")) +
  geom_ribbon(aes(x = NO3adj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = NO3adj, y = fit), inherit.aes=FALSE) 
EffNO3

EffNH4 <- ggplot(testCINH4, 
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
        legend.position = "none",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIP ("*mu*"M)")) 
EffNH4

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
  labs(color = expression("Species")) 
EffP
grid.arrange(EffNO3,EffNH4,EffP,ncol=3)
