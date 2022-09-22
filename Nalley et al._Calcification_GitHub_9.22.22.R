#### NALLEY ET AL. 2022 NUTRIENT STRESSOR THRESHOLDS ####
### Calcification (mg CaCO3 cm-2 day-1) ###

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
calc <- as.data.frame(read.csv("ALL-Calcification.csv", header=TRUE))
str(calc) ## everything is factor vs. numeric

# creating new variable that has both Genus and species
calc <- calc %>%
  mutate(Gsp = paste(Current.genus.name, Current.species.name, sep = " ")) 
head(calc)
### estimating ambient concentrations of N & P if listed as zero
calc$Nadj <- ifelse(calc$TOTAL.N == 0, 0.1, calc$TOTAL.N)
calc$Padj <- ifelse(calc$PO4 == 0, 0.02, calc$PO4)
### adding log variables
calc$log10N <- log10(calc$Nadj+1)
calc$log10P <- log10(calc$Padj+1)

### ordering data before calculating hedges
str(calc)
calc.ord <- calc[order(calc$Reference, calc$Control),]

### calculating hedges D 
covar <- by(calc.ord, calc.ord$Reference, function(x) 
  covar.smd(Response.level, SD, N, "smd", method="hedges", data = x))
calc.ord$smd <- unlist(lapply(covar, function(x) x$y))
calc.ord$vmd <- unlist(lapply(covar, function(x) x$v))
str(calc.ord) 
test <- calc.ord[,c("Reference", "Control", "smd")]

### removing controls for models
calc2 <- subset(calc.ord, Control=="exp")
summary(calc2)
calc2$smdpos <- as.factor(ifelse(calc2$smd > 0, "Increase", "Decrease"))

str(calc2)
## 21 Studies with 44 References within
CalcTreatSummary <- calc2 %>% 
  summarize(Nmin = min(Nadj),
            Nmax = max(Nadj),
            Pmin = min(Padj),
            Pmax = max(Padj),
            DurMean = mean(Exposure.in.days),
            DurMin = min(Exposure.in.days),
            DurMax = max(Exposure.in.days))

## dotplot
ggplot(calc2, 
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
  scale_size(range = c(3,15)) +
  geom_jitter(alpha = .5, position = position_jitter(width = .07, height = 0.03)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = expression("DIP ("*mu*")M"),
       size = "Treatment\nEffect Size", color = "Calc. Rate") +
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
                data = calc2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod1) 

## log10N + log10P with simple random
mod2 <- mixmeta(smd ~ log10N + log10P, 
                random = ~ 1 | Reference, 
                data = calc2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod2) 

AIC(mod1, mod2)
BIC(mod1, mod2)

## log10N with poly2 + log10P with poly2
mod3 <- mixmeta(smd ~ poly(log10N, 2, raw = TRUE) + poly(log10P, 2, raw = TRUE), 
                random = ~ 1 | Reference, 
                data = calc2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod3)

## log10N with poly2 + log10P
mod4 <- mixmeta(smd ~ poly(log10N, 2, raw = TRUE) + log10P, 
                random = ~ 1 | Reference, 
                data = calc2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod4) 

## log10N with poly3+ log10P with poly3
mod5 <- mixmeta(smd ~ poly(log10N, 3, raw = TRUE) + poly(log10P, 3, raw = TRUE),
                random = ~ 1 | Reference, 
                data = calc2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod5)

## log10N with poly3 + log10P 
mod6 <- mixmeta(smd ~ poly(log10N, 3, raw = TRUE) + log10P, 
                random = ~ 1 | Reference, 
                data = calc2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod6) 

AIC(mod2, mod3, mod4, mod5, mod6) 
BIC(mod2, mod3, mod4, mod5, mod6) 

## log10N with poly3 + log10P + Exposure in days
mod7 <- mixmeta(smd ~ log10N + log10P + Gsp, 
                random = ~ 1 | Reference, 
                data = calc2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod7) 

mod8 <- mixmeta(smd ~ log10N + log10P + Exposure.in.days, 
                 random = ~ 1 | Reference, 
                 data = calc2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod8) 
AIC(mod2, mod7, mod8) 
BIC(mod2, mod7, mod8) 


### PLOTTING: PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS ####
newdata.predictN <- data.frame(log10N = calc2$log10N,
                               log10P = median(calc2$log10P))
newdata.predictP <- data.frame(log10N = median(calc2$log10N),
                               log10P = calc2$log10P)
predN <- predict(mod2, newdata=newdata.predictN, ci=TRUE)
predP <- predict(mod2, newdata=newdata.predictP, ci=TRUE)
testCIN <- cbind(calc2, predN)
testCIP <- cbind (calc2, predP)

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
        legend.position = c(0.15, 0.2),
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
        legend.position = c(0.85, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIN ("*mu*"M)")) 
EffP
grid.arrange(EffN,EffP,ncol=2)
summary(mod2)


#### and colored by species for supplemental ####
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
        legend.position = c(0.85, 0.25),
        legend.text=element_text(size=6, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("")) 
EffP
grid.arrange(EffN,EffP,ncol=2)


#### and colored by family for supplemental ####
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
        legend.position = c(0.85, 0.25),
        legend.text=element_text(size=8),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("")) 
EffP
grid.arrange(EffN,EffP,ncol=2)

#### and colored by morphology for supplemental ####
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
        legend.position = c(0.85, 0.25),
        legend.text=element_text(size=8),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("")) 
EffP
grid.arrange(EffN,EffP,ncol=2)

