#### NALLEY ET AL. 2022 NUTRIENT THRESHOLD ANALYSIS #####
### Growth ###

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
library(ggpubr)
library(ggplot2)
library(gridExtra)

breaks_log10 <- function(x) {
  low <- floor(log10(min(x)))
  high <- ceiling(log10(max(x)))
  
  10^(seq.int(low, high))
}

getwd()
setwd("/")
grow <- as.data.frame(read.csv("ALL-Growth.csv", header=TRUE))
sub.grow <- as.data.frame(read.csv("ALL-Growth_sub.csv", header=TRUE))
str(grow) ## everything is factor vs. numeric
str(sub.grow)
summary(sub.grow)

# creating new variable that has both Genus and species
grow <- grow %>%
  mutate(Gsp = paste(Genus, Species, sep = " ")) 
head(grow)
sub.grow <- sub.grow %>%
  mutate(Gsp = paste(Genus, Species, sep = " ")) 

### estimating ambient concentrations of N & P if listed as zero
grow$Nadj <- ifelse(grow$TOTAL.N == 0, 0.1, grow$TOTAL.N)
grow$Padj <- ifelse(grow$PO4 == 0, 0.02, grow$PO4)
sub.grow$Nadj <- ifelse(sub.grow$TOTAL.N == 0, 0.1, sub.grow$TOTAL.N)
sub.grow$Padj <- ifelse(sub.grow$PO4 == 0, 0.02, sub.grow$PO4)
sub.grow$NO3adj <- ifelse(sub.grow$NO3 == 0, 0.1, sub.grow$NO3)
sub.grow$NH4adj <- ifelse(sub.grow$NH4 == 0, 0.1, sub.grow$NH4)
### adding log variables
grow$log10N <- log10(grow$Nadj+1)
grow$log10P <- log10(grow$Padj+1)
sub.grow$log10N <- log10(sub.grow$Nadj+1)
sub.grow$log10P <- log10(sub.grow$Padj+1)
sub.grow$log10NO3 <- log10(sub.grow$NO3adj+1)
sub.grow$log10NH4 <- log10(sub.grow$NH4adj+1)

### ordering data before calculating hedges
str(grow)
grow.ord <- grow[order(grow$Reference, grow$Control),]
sub.grow.ord <- sub.grow[order(sub.grow$Reference, sub.grow$Control),]

### calculating hedges D 
covar <- by(grow.ord, grow.ord$Reference, function(x) 
  covar.smd(Response.level, SD, N, "smd", method="hedges", data = x))
grow.ord$smd <- unlist(lapply(covar, function(x) x$y))
grow.ord$vmd <- unlist(lapply(covar, function(x) x$v))
str(grow.ord) 
test <- grow.ord[,c("Reference", "Control", "smd")]
### removing controls for models
grow2 <- subset(grow.ord, Control=="exp")
summary(grow2)
## & adding in smdpos factor
grow2$smdpos <- as.factor(ifelse(grow2$smd > 0, "Increase", "Decrease"))
str(grow2)

## & without Rasmussen 1994
### calculating hedges D 
sub.covar <- by(sub.grow.ord, sub.grow.ord$Reference, function(x) 
  covar.smd(Response.level, SD, N, "smd", method="hedges", data = x))
sub.grow.ord$smd <- unlist(lapply(sub.covar, function(x) x$y))
sub.grow.ord$vmd <- unlist(lapply(sub.covar, function(x) x$v))
sub.grow.ord <- as.data.frame(sub.grow.ord)
str(sub.grow.ord) 
test <- sub.grow.ord[,c("Reference", "Control", "smd")]
### removing controls for models
sub.grow2 <- subset(sub.grow.ord, Control=="exp")
summary(sub.grow2)
## & adding in smdpos factor
sub.grow2$smdpos <- as.factor(ifelse(sub.grow2$smd > 0, "Increase", "Decrease"))
str(sub.grow2)

## 7 Studies with 9 References within
GrowTreatSummary <- grow2 %>% 
  summarize(Nmin = min(Nadj),
            Nmax = max(Nadj),
            Pmin = min(Padj),
            Pmax = max(Padj),
            DurMean = mean(Duration.in.days),
            DurMin = min(Duration.in.days),
            DurMax = max(Duration.in.days))

## dotplot
ggplot(sub.grow2, 
       aes(x = Nadj, y = Padj, 
           color = smdpos,
           size = abs(smd))) + #geom_pointrange() +
  ## surface water from station aloha 2019 annual report (surface ocean)
  annotate("point", x = 0.03, y = 0.03, colour = "gray", size = 4, shape = 8) +
  ## ambient reef from Nakajima et al. 2015 (Malaysia)
  annotate("point", x = 0.75, y = 0.1, colour = "gray",  size = 4, shape = 8) +
  ## ambient reef from Silbiger et al. 2018 (Kaneohe Bay)
  annotate("point", x = .15, y = .15, colour = "gray",  size = 4, shape = 8) + 
  ## ambient (but high) from AIMS reef in Fabricius et al. 2013 (GBR)
  annotate("point", x = .44, y = .38, colour = "gray",  size = 4, shape = 8) + 
  ## high reef from Silbiger et al. 2018 (throughout Pacific)
  annotate("point", x = 7.6, y = 2.6, colour = "gray",  alpha = 0.8, size = 4, shape = 8) + 
  ## highest measured at SGD seep in Lubarsky et al. 2018 (Maunalua Bay)
  annotate("point", x = 32.39, y = 0.89, colour = "gray",  size = 4, shape = 8) + 
  scale_x_log10(breaks = breaks_log10) +  
  scale_y_log10(breaks = breaks_log10) + annotation_logticks() +
  scale_size(range = c(3,15)) +
  geom_jitter(alpha = .5, position = position_jitter(width = .07, height = 0.03)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20)) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = expression("DIP ("*mu*"M)"),
       size = "Treatment\nEffect Size", color = "Growth") + 
  guides(color = guide_legend(order = 2, override.aes = list(size=4)),
         size = guide_legend(order = 1)) 


### COVARIANCE MATRICES
newlist <- list(NA)
for (i in seq(1,length(covar))) {
  newlist[i] <- list(covar[[i]]$S)
}
##
sub.newlist <- list(NA)
for (i in seq(1,length(sub.covar))) {
  sub.newlist[i] <- list(sub.covar[[i]]$S)
}

### LOOKING AT SMD
## N + P with simple random
mod1 <- mixmeta(smd ~ Nadj + Padj, 
                random = ~ 1 | Reference, 
                data = sub.grow2, method = "ml", 
                control = list(addSlist = sub.newlist))
summary(mod1)

## log10N + log10P with simple random
mod2 <- mixmeta(smd ~ log10N + log10P, 
                random = ~ 1 | Reference, 
                data = sub.grow2, method = "ml", 
                control = list(addSlist = sub.newlist))
summary(mod2) 

AIC(mod1, mod2)
BIC(mod1, mod2)

## Nadj with poly2 + Padj with poly2
mod3 <- mixmeta(smd ~ poly(Nadj, 2, raw = TRUE) + poly(Padj, 2, raw = TRUE), 
                random = ~ 1 | Reference, 
                data = sub.grow2, method = "ml", 
                control = list(addSlist = sub.newlist))
summary(mod3) # 

## Padj with poly2 + Nadj
mod4 <- mixmeta(smd ~ poly(Padj, 2, raw = TRUE) + Nadj, 
                random = ~ 1 | Reference, 
                data = sub.grow2, method = "ml", 
                control = list(addSlist = sub.newlist))
summary(mod4) # 

## Nadj with poly3+ Padj with poly3
mod5 <- mixmeta(smd ~ poly(Nadj, 3, raw = TRUE) + poly(Padj, 3, raw = TRUE),
                random = ~ 1 | Reference, 
                data = sub.grow2, method = "ml", 
                control = list(addSlist = sub.newlist))
summary(mod5) # 

## Padj with poly3 + Nadj 
mod6 <- mixmeta(smd ~ poly(Padj, 3, raw = TRUE) + Nadj, 
                random = ~ 1 | Reference, 
                data = sub.grow2, method = "ml", 
                control = list(addSlist = sub.newlist))
summary(mod6) # 

AIC(mod1, mod3, mod4, mod5, mod6) ##
BIC(mod1, mod3, mod4, mod5, mod6) ## 

## Nadj + Padj + Gsp
mod7 <- mixmeta(smd ~ Nadj + Padj + Gsp, 
                random = ~ 1 | Reference, 
                data = sub.grow2, method = "ml", 
                control = list(addSlist = sub.newlist))
summary(mod7) ##

## Nadj + Padj + Exp Duration
mod8 <- mixmeta(smd ~ Nadj + Padj + Duration.in.days, 
                random = ~ 1 | Reference, 
                data = sub.grow2, method = "ml", 
                control = list(addSlist = sub.newlist))
summary(mod8) ## 

AIC(mod1, mod7, mod8)
BIC(mod1, mod7, mod8)


### PLOTTING: PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS ####
newdata.predictN <- data.frame(Nadj = sub.grow2$Nadj,
                               Padj = median(sub.grow2$Padj),
                               Duration.in.days = median(sub.grow2$Duration.in.days))
newdata.predictP <- data.frame(Nadj = median(sub.grow2$Nadj),
                               Padj = sub.grow2$Padj,
                               Duration.in.days = median(sub.grow2$Duration.in.days))
newdata.predictE <- data.frame(Nadj = median(sub.grow2$Nadj),
                               Padj = median(sub.grow2$Padj),
                               Duration.in.days = sub.grow2$Duration.in.days)

predN <- predict(mod8, newdata=newdata.predictN, ci=TRUE)
predP <- predict(mod8, newdata=newdata.predictP, ci=TRUE)
predE <- predict(mod8, newdata=newdata.predictE, ci=TRUE)
testCIN <- cbind(sub.grow2, predN)
testCIP <- cbind (sub.grow2, predP)
testCIE <- cbind (sub.grow2, predE)

EffN <- ggplot(testCIN, 
               aes(x =Nadj,
                   y = smd, color = Padj,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) + ylim(-1.2,3.4) +
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
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = Nadj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP, 
               aes(x =Padj,
                   y = smd, color = Nadj,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) + ylim(-1.2,3.4) +
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

## exposure duration as a predictor
EffE <- ggplot(testCIE, 
               aes(x = Duration.in.days,
                   y = smd, color = Padj,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) + ylim(-1.2,3.4) +
  labs(x = expression("Exposure (in days)"),
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
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = Duration.in.days, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Duration.in.days, y = fit), inherit.aes=FALSE) 
EffE
grid.arrange(EffN,EffP,EffE,ncol=3)
summary(mod8)

#### colored by Gsp for supplemental ####
EffN <- ggplot(testCIN, 
               aes(x =Nadj,
                   y = smd, color = Gsp,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) + ylim(-1.2,3.4) +
  labs(x = expression("DIN ("*mu*"M)"),
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
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = Nadj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP, 
               aes(x =Padj,
                   y = smd, color = Gsp,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) + ylim(-1.2,3.4) +
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
        legend.position = c(0.25, 0.85),
        legend.text=element_text(size=10, face = "italic"),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("")) +
  geom_ribbon(aes(x = Padj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Padj, y = fit), inherit.aes=FALSE) 
EffP

## exposure duration as a predictor
EffE <- ggplot(testCIE, 
               aes(x = Duration.in.days,
                   y = smd, color = Gsp,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) + ylim(-1.2,3.4) +
  labs(x = expression("Exposure (in days)"),
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
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = Duration.in.days, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Duration.in.days, y = fit), inherit.aes=FALSE) 
EffE
grid.arrange(EffN,EffP,EffE,ncol=3)

#### colored by Family for supplemental ####
EffN <- ggplot(testCIN, 
               aes(x =Nadj,
                   y = smd, color = Family,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) + ylim(-1.2,3.4) +
  labs(x = expression("DIN ("*mu*"M)"),
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
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = Nadj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP, 
               aes(x =Padj,
                   y = smd, color = Family,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) + ylim(-1.2,3.4) +
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
        legend.position = c(0.25, 0.85),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("")) +
  geom_ribbon(aes(x = Padj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Padj, y = fit), inherit.aes=FALSE) 
EffP

## exposure duration as a predictor
EffE <- ggplot(testCIE, 
               aes(x = Duration.in.days,
                   y = smd, color = Family,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) + ylim(-1.2,3.4) +
  labs(x = expression("Exposure (in days)"),
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
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = Duration.in.days, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Duration.in.days, y = fit), inherit.aes=FALSE) 
EffE
grid.arrange(EffN,EffP,EffE,ncol=3)



#### colored by morphology for supplemental ####
EffN <- ggplot(testCIN, 
               aes(x =Nadj,
                   y = smd, color = Coral.Type,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) + ylim(-1.2,3.4) +
  labs(x = expression("DIN ("*mu*"M)"),
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
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = Nadj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) 
EffN

## and P as a predictor
EffP <- ggplot(testCIP, 
               aes(x =Padj,
                   y = smd, color = Coral.Type,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) + ylim(-1.2,3.4) +
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
        legend.position = c(0.15, 0.9),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("")) +
  geom_ribbon(aes(x = Padj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Padj, y = fit), inherit.aes=FALSE) 
EffP

## exposure duration as a predictor
EffE <- ggplot(testCIE, 
               aes(x = Duration.in.days,
                   y = smd, color = Coral.Type,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) + ylim(-1.2,3.4) +
  labs(x = expression("Exposure (in days)"),
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
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = Duration.in.days, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Duration.in.days, y = fit), inherit.aes=FALSE) 
EffE
grid.arrange(EffN,EffP,EffE,ncol=3)

