#### Nutrients Analysis #####
### All-Calc ###
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
setwd("/Users/eileennalley/Desktop/WORK/Thresholds/Nutrients/Analysis")
calc <- as.data.frame(read.csv("ALL-Calcification.csv", header=TRUE))
str(calc) ## everything is factor vs. numeric

# creating new variable that has both Genus and species
calc <- calc %>%
  mutate(Gsp = paste(Current.genus.name, Current.species.name, sep = "_")) 
head(calc)
### estimating ambient concentrations of N & P if listed as zero
calc$Nadj <- ifelse(calc$TOTAL.N == 0, 0.1, calc$TOTAL.N)
calc$Padj <- ifelse(calc$PO4 == 0, 0.02, calc$PO4)
### adding log variables
calc$log10N <- log10(calc$Nadj)
calc$log10P <- log10(calc$Padj)

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


CalcN <- ggplot(calc, 
       aes(x = Nadj, y = Response.level, 
           ymin = Response.level-SE,
           ymax = Response.level+SE,
           color = Control)) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +  
  geom_smooth(method = "loess", color = "darkgrey", fill = "lightgrey") +
  geom_pointrange(alpha = 0.7) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.75, 0.85),
        legend.background = element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=12)) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = expression("Calcification (mg "~CaCO[3]*""~cm^-2*""~day^-1*")"),
       color = element_blank()) +
  scale_color_brewer(palette="Set2")
CalcN

CalcP <- ggplot(calc, 
       aes(x = Padj, y = Response.level, 
           ymin = Response.level-SE,
           ymax = Response.level+SE,
           color = Control)) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") + 
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
       y = element_blank(),
       color = element_blank()) +
  scale_color_brewer(palette="Set2")
CalcP
grid.arrange(CalcN,CalcP,ncol=2)



calc3 <- calc2%>%mutate(NadjQuantBins = cut(Nadj, 
                                            breaks = unique(quantile(Nadj,probs=seq.int(0,1, by=1/4))), 
                                            include.lowest=TRUE))
calc3 <- calc3%>%mutate(PadjQuantBins = cut(Padj, 
                                            breaks = unique(quantile(Padj,probs=seq.int(0,1, by=1/4))), 
                                            include.lowest=TRUE))
N <- ggplot(calc3, aes(x = NadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Proportion of Reported Change in Calcification",
       fill = "") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(size = 15)) +
  scale_x_discrete(labels = c("[0.2,1.56]" = "0.2-1.6",
                              "(1.56,5.2]" = "1.6-5.2",
                              "(5.2,20]" = "5.2-20",
                              "(20,50]" = "20-50")) 
P <- ggplot(calc3, aes(x = PadjQuantBins, fill = smdpos)) +
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
        axis.text.x = element_text(size = 15)) +
  scale_x_discrete(labels = c("[0.02,0.04]" = "0.02-0.04",
                              "(0.04,0.12]" = "0.04-0.1",
                              "(0.12,5]" = "0.1-5.0")) 
grid.arrange(N,P,ncol=2)


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

plot(smd~log10N, data = calc2)
plot(smd~log10P, data = calc2)
### LOOKING AT SMD
## N + P with zero
mod5 <- mixmeta(smd ~ 0 + Nadj + Padj, 
                random = ~ 1 | Reference, 
                data = calc2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod5) # I2 = 61.1

## N + P + N*P with zero
mod6 <- mixmeta(smd ~ 0 + Nadj + Padj + Nadj:Padj, 
                random = ~ 1| Reference, 
                data = calc2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod6) # I2 = 61.3

## log10N + log10P with zero
mod7 <- mixmeta(smd ~ 0 + log10N + log10P, 
                random = ~ 1 | Reference, 
                data = calc2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod7) # I2 = 59.9

## log10N + log10P + log10N*log10P with zero
mod8 <- mixmeta(smd ~ 0 + log10N + log10P + log10N:log10P, 
                random = ~ 1 | Reference, 
                data = calc2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod8) # I2 = 60.9

## N + P with simple random
mod9 <- mixmeta(smd ~ Nadj + Padj, 
                random = ~ 1 | Reference, 
                data = calc2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod9) # I2 = 56.3
mod9b <- mixmeta(smd ~ Nadj + Padj + Exposure.in.days, 
                random = ~ 1 | Reference, 
                data = calc2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod9b) #
AIC(mod9, mod9b)

## N + P + N*P with simple random
mod10 <- mixmeta(smd ~ Nadj + Padj + Nadj:Padj, 
                 random = ~ 1| Reference, 
                 data = calc2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod10) # I2 = 57.3

## log10N + log10P with simple random
mod11 <- mixmeta(smd ~ log10N + log10P, 
                 random = ~ 1 | Reference, 
                 data = calc2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod11) # I2 = 55.9

mod11b <- mixmeta(smd ~ log10N + log10P + Exposure.in.days, 
                 random = ~ 1 | Reference, 
                 data = calc2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod11b)

mod11c <- mixmeta(smd ~ log10N + log10P + Gsp, 
                  random = ~ 1 | Reference, 
                  data = calc2, method = "ml", 
                  control = list(addSlist = newlist))
summary(mod11c)

mod11d <- mixmeta(smd ~ log10N + log10P + Current.genus.name, 
                  random = ~ 1 | Reference, 
                  data = calc2, method = "ml", 
                  control = list(addSlist = newlist))
summary(mod11d)

mod11e <- mixmeta(smd ~ log10N + log10P + Author.s..Year, 
                  random = ~ 1 | Reference, 
                  data = calc2, method = "ml", 
                  control = list(addSlist = newlist))
summary(mod11e)
AIC(mod11, mod11b, mod11c, mod11d, mod11e)

## log10N + log10 P + log10N*log10P with simple random
mod12 <- mixmeta(smd ~ log10N + log10P + log10N:log10P, 
                 random = ~ 1 | Reference, 
                 data = calc2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod12) # I2 = 56.4
AIC(mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12) # mod11,9
BIC(mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12) # mod11,9
summary(mod11)
resid <- resid(mod11)
fitted <- fitted(mod11)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)
## pretty good residuals 

## generating 1st & 3rd Q for nonlinear models
summary(calc2$log10N)
summary(calc2$log10P)
summary(calc2$Nadj)
summary(calc2$Padj)

## nonlin
mod13 <- mixmeta(smd ~ 0 + ns(log10N, knots = c(0.1928, 1.3010)) +
                   ns(log10P, knots = c(-1.6990, -0.9208)),
                 random = ~ 1 | Reference, 
                 data = calc2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod13) # I2 = 57.6

mod13b <- mixmeta(smd ~ 0 + ns(Nadj, knots = c(1.56, 20.00)) +
                    ns(Padj, knots = c(0.0200, 0.1200)),
                  random = ~ 1 | Reference, 
                  data = calc2, method = "ml", 
                  control = list(addSlist = newlist))
summary(mod13b) # I2 = 59

## with quartiles
mod14quart <- mixmeta(smd ~ 0+ ns(log10N, knots=c(0.1928,0.7160,1.3010)) + 
                        ns(log10P, knots=c(-1.6990,-1.3979,-0.9208)), 
                      data=calc2,
                      random= ~ 1|Reference, 
                      method="ml", 
                      control=list(addSlist=newlist))
summary(mod14quart) # I2 = 52

mod14quartb <- mixmeta(smd ~ 0+ ns(Nadj, knots=c(1.56,5.20,20)) + 
                         ns(Padj, knots=c(.02,0.0400,0.1200)), 
                       data=calc2,
                       random= ~ 1|Reference, 
                       method="ml", 
                       control=list(addSlist=newlist))
summary(mod14quartb) # I2 = 57.1
AIC(mod11, mod13, mod13b, mod14quart, mod14quartb)
## mod11 best 

summary(mod11)
# PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS
newdata.predictN <- data.frame(log10N = calc2$log10N,
                               log10P = median(calc2$log10P))
newdata.predictP <- data.frame(log10N = median(calc2$log10N),
                               log10P = calc2$log10P)
predN <- predict(mod11, newdata=newdata.predictN, ci=TRUE)
predP <- predict(mod11, newdata=newdata.predictP, ci=TRUE)
testCIN <- cbind(calc2, predN)
testCIP <- cbind (calc2, predP)
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
        legend.position = c(0.85, 0.85),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
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
        legend.position = c(0.85, 0.85),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIN ("*mu*"M)")) 

grid.arrange(ESNlin,ESPlin,ncol=2)
summary(mod11)

