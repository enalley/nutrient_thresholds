#### NALLEY ET AL. 2022 NUTRIENT STRESSOR THRESHOLDS ####
### Growth (mm day-1) ###

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
grow <- as.data.frame(read.csv("ALL-Growth_sub.csv", header=TRUE))
str(grow) ## everything is factor vs. numeric

# creating new variable that has both Genus and species
grow <- grow %>%
  mutate(Gsp = paste(Genus, Species, sep = "_")) 
head(grow)

### estimating ambient concentrations of N & P if listed as zero
grow$Nadj <- ifelse(grow$TOTAL.N == 0, 0.1, grow$TOTAL.N)
grow$Padj <- ifelse(grow$PO4 == 0, 0.02, grow$PO4)
### adding log variables
grow$log10N <- log10(grow$Nadj)
grow$log10P <- log10(grow$Padj)

### ordering data before calculating hedges
grow.ord <- grow[order(grow$Reference, grow$Control),]


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

## 7 Studies with 9 References within
GrowTreatSummary <- grow2 %>% 
  summarize(Nmin = min(Nadj),
            Nmax = max(Nadj),
            Pmin = min(Padj),
            Pmax = max(Padj),
            DurMean = mean(Duration.in.days),
            DurMin = min(Duration.in.days),
            DurMax = max(Duration.in.days))


GrowN <- ggplot(grow, 
       aes(x = Nadj, y = Response.level, 
           ymin = Response.level-SE,
           ymax = Response.level+SE,
           color = Control)) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") + ylim(-0.02,0.25) +
  geom_smooth(method = "loess", color = "darkgrey", fill = "lightgrey") +
  geom_pointrange(alpha = 0.7) +  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.22, 0.83),
        legend.background = element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=12)) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Growth (mm "~day^-1*")",
       color = element_blank()) +
  scale_color_brewer(palette="Set2")
GrowN

GrowP <- ggplot(grow, 
                aes(x = Padj, y = Response.level, 
                    ymin = Response.level-SE,
                    ymax = Response.level+SE,
                    color = Control)) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") + ylim(-.02,0.25) +
  geom_smooth(method = "loess", color = "darkgrey", fill = "lightgrey") +
  geom_pointrange(alpha = 0.7) + theme_bw() + 
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
       y = "Growth (mm "~day^-1*")",
       color = element_blank()) +
  scale_color_brewer(palette="Set2")
GrowP
grid.arrange(GrowN,GrowP,ncol=2)


grow3 <- grow2%>%mutate(NadjQuantBins = cut(Nadj, 
                                                breaks = unique(quantile(Nadj,probs=seq.int(0,1, by=1/4))), 
                                                include.lowest=TRUE))
grow3 <- grow3%>%mutate(PadjQuantBins = cut(Padj, 
                                                breaks = unique(quantile(Padj,probs=seq.int(0,1, by=1/4))), 
                                                include.lowest=TRUE))
N <- ggplot(grow3, aes(x = NadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("Total N "*mu*"M"),
       y = "Proportion of Reported Change in Growth",
       fill = "") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(size = 15)) +
  scale_x_discrete(labels = c("[0.1,2.23]" = "0.1-2.2",
                              "(2.23,15.1]" = "2.2-15.1",
                              "(15.1,20.2]" = "15.1-20.2",
                              "(20.2,50]" = "20.2-50"))
P <- ggplot(grow3, aes(x = PadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression(""~PO[4]*" "*mu*"M"),
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
        axis.text.x = element_text(size = 15)) +
  scale_x_discrete(labels = c("[0.02,0.39]" = "0.02-0.4",
                              "(0.39,4.85]" = "0.4-4.9",
                              "(4.85,16.1]" = "4.9-16.1"))
grid.arrange(N,P,ncol=2)



## dotplot
ggplot(grow2, 
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


grow3 <- grow2%>%mutate(NadjQuantBins = cut(Nadj, 
                                          breaks = unique(quantile(Nadj,probs=seq.int(0,1, by=1/4))), 
                                          include.lowest=TRUE))
grow3 <- grow3%>%mutate(PadjQuantBins = cut(Padj, 
                                          breaks = unique(quantile(Padj,probs=seq.int(0,1, by=1/4))), 
                                          include.lowest=TRUE))

N <- ggplot(grow3, aes(x = NadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Proportion of Reported Change in Growth",
       fill = "Growth Rate") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none") +
  scale_x_discrete(labels = c("[0.1,4.58]" = "0.1-4.6",
                              "(4.58,20]" = "4.6-20",
                              "(20,28.9]" = "20-28.9",
                              "(28.9,50]" = "28.9-50"))

P <- ggplot(grow3, aes(x = PadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Proportion",
       fill = "Growth Rate") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.position = "none")+
  scale_x_discrete(labels = c("[0.02,0.16]" = "0.02-0.2",
                              "(0.16,3.57]" = "0.2-3.6",
                              "(3.57,16.1]" = "3.6-16.1"))
grid.arrange(N,P,ncol=2)

### COVARIANCE MATRICES
newlist <- list(NA)
for (i in seq(1,length(covar))) {
  newlist[i] <- list(covar[[i]]$S)
}


### LOOKING AT SMD
## N + P with zero
mod5 <- mixmeta(smd ~ 0 + Nadj + Padj, 
                random = ~ 1 | Reference, 
                data = grow2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod5) # I2 = 75.5
mod5b <- mixmeta(smd ~ 0 + Nadj + Padj + Duration.in.days, 
                random = ~ 1 | Reference, 
                data = grow2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod5b)
AIC(mod5, mod5b)
mod5c <- mixmeta(smd ~ 0 + Nadj + Padj + I(Padj^3), 
                 random = ~ 1 | Reference, 
                 data = grow2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod5c)

## N + P + N*P with zero
mod6 <- mixmeta(smd ~ 0 + Nadj + Padj + Nadj:Padj, 
                random = ~ 1| Reference, 
                data = grow2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod6) # I2 = 76.5

## log10N + log10P with zero
mod7 <- mixmeta(smd ~ 0 + log10N + log10P, 
                random = ~ 1 | Reference, 
                data = grow2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod7) # I2 = 79.7

## log10N + log10P + log10N*log10P with zero
mod8 <- mixmeta(smd ~ 0 + log10N + log10P + log10N:log10P, 
                random = ~ 1 | Reference, 
                data = grow2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod8) # I2 = 80.4

## N + P with simple random
mod9 <- mixmeta(smd ~ Nadj + Padj, 
                random = ~ 1 | Reference, 
                data = grow2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod9) # I2 = 76.2

## N + P + N*P with simple random
mod10 <- mixmeta(smd ~ Nadj + Padj + Nadj:Padj, 
                 random = ~ 1| Reference, 
                 data = grow2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod10) # I2 = 77.3

## log10N + log10P with simple random
mod11 <- mixmeta(smd ~ log10N + log10P, 
                 random = ~ 1 | Reference, 
                 data = grow2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod11) # I2 = 78.9

## log10N + log10 P + log10N*log10P with simple random
mod12 <- mixmeta(smd ~ log10N + log10P + log10N:log10P, 
                 random = ~ 1 | Reference, 
                 data = grow2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod12) # I2 = 79.6
AIC(mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12) # mod5
BIC(mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12) # mod5
summary(mod5)
resid <- resid(mod5)
fitted <- fitted(mod5)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)

## generating 1st & 3rd Q for nonlinear models
summary(grow2$log10N)
summary(grow2$log10P)
summary(grow2$Nadj)
summary(grow2$Padj)

## nonlin
## simple random
mod14 <- mixmeta(smd ~ 0 + ns(log10N, knots=c(0.6577,1.4400)) + 
                   ns(log10P, knots=c(-1.6990,0.5060)), 
                 data=grow2,
                 random= ~ 1| Reference, 
                 method="ml", 
                 control=list(addSlist=newlist))
summary(mod14) # I2 = 70.6
mod14b <- mixmeta(smd ~ 0 + ns(Nadj, knots = c(4.576, 28.875)) +
                    ns(Padj, knots = c(0.0200, 3.570)), 
                  data=grow2,
                  random= ~ 1| Reference, 
                  method="ml", 
                  control=list(addSlist=newlist))
summary(mod14b)

## with quartiles
mod14quart <- mixmeta(smd ~ 0+ ns(log10N, knots=c(0.6577,1.3010,1.4400)) + 
                        ns(log10P, knots=c(-1.6990,-0.7959,-0.5060)), 
                      data=grow2,
                      random= ~ 1|Reference, 
                      method="ml", 
                      control=list(addSlist=newlist))
summary(mod14quart) # I2 = 71.3
mod14quartb <- mixmeta(smd ~ 0+ ns(Nadj, knots=c(4.576,20,28.875)) + 
                         ns(Padj, knots=c(.02,0.160,3.570)), 
                       data=grow2,
                       random= ~ 1|Reference, 
                       method="ml", 
                       control=list(addSlist=newlist))
summary(mod14quartb)

AIC(mod5, mod9, mod14, mod14b, mod14quart, mod14quartb)
BIC(mod5, mod9, mod14, mod14b, mod14quart, mod14quartb)
## mod 9 is the best but the relationship with PO4 is clearly nonlinear so using mod14b which is the next best
summary(mod14b)
resid <- resid(mod14b)
fitted <- fitted(mod14b)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)

# PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS
## with linear model 3
newdata.predictN <- data.frame(Nadj = grow2$Nadj,
                               Padj = median(grow2$Padj))
newdata.predictP <- data.frame(Nadj = median(grow2$Nadj),
                               Padj = grow2$Padj)
predN <- predict(mod14b, newdata=newdata.predictN, ci=TRUE)
predP <- predict(mod14b, newdata=newdata.predictP, ci=TRUE)
testCIN <- cbind(grow2, predN)
testCIP <- cbind (grow2, predP)
min_conc_N <- testCIN %>% 
  mutate(overlap0 = 0 >= ci.lb & 0 <= ci.ub) %>% 
  filter(overlap0==FALSE) %>% 
  summarize(min_N = min(Nadj),
            max_N = max(Nadj),
            min_P = min(Padj),
            max_P = max(Padj))

## and now plotting for nitrogen as a predictor
ESNlin <- ggplot(testCIN, 
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
  annotation_logticks(sides = "b") + ylim(-1.5,4) +
  labs(color = expression("DIP ("*mu*"M)")) 
ESNlin

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
        legend.position = c(0.15, 0.85),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") + ylim(-1.5,4) +
  labs(color = expression("DIN ("*mu*"M)")) +
  geom_ribbon(aes(x = Padj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Padj, y = fit), inherit.aes=FALSE) 
grid.arrange(ESNlin,ESPlin,ncol=2)

