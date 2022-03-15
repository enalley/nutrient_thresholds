#### NALLEY ET AL. 2022 NUTRIENT THRESHOLD ANALYSIS #### 
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

symb <- as.data.frame(read.csv("ALL-CellsSA.csv", header=TRUE))
str(symb) ## everything is correct as factor vs. numeric
symb <- symb[,c(1:17)]
symb <- as.data.frame(symb[!symb$Reference == "NU20",])

str(symb)
## only press experiments

# creating new variable that has both Genus and species
symb <- symb %>%
  mutate(Gsp = paste(Genus, Species, sep = "_")) 
head(symb)
### estimating ambient concentrations of N & P if listed as zero
symb$Nadj <- ifelse(symb$TOTAL.N == 0, 0.1, symb$TOTAL.N)
symb$Padj <- ifelse(symb$PO4 == 0, 0.02, symb$PO4)
### adding log variables
symb$log10N <- log10(symb$Nadj)
symb$log10P <- log10(symb$Padj)

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


SymbN <- ggplot(symb, 
                aes(x = Nadj, y = Response.level, 
                    ymin = Response.level-SE,
                    ymax = Response.level+SE,
                    color = Control)) + 
  geom_hline(yintercept = 3, color = "lightgrey", linetype = "dashed")+
  geom_smooth(method = "loess", color = "darkgrey", fill = "lightgrey") +
  geom_pointrange(alpha = 0.7) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +  ylim(-0.5,5.25) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.15, 0.9),
        legend.background = element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=12)) +
  labs(x = expression("DIN ("*mu*"M)"),
       y = expression("Symbiont Density ("~10^6*" cells"~cm^-2*")"),
       color = element_blank()) +
  scale_color_brewer(palette="Set2")
SymbN

SymbP <- ggplot(symb, 
                aes(x = Padj, y = Response.level, 
                    ymin = Response.level-SE,
                    ymax = Response.level+SE,
                    color = Control)) +
  geom_hline(yintercept = 3, color = "lightgrey", linetype = "dashed")+
  geom_smooth(method = "loess", color = "darkgrey", fill = "lightgrey") +
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +  ylim(-0.5,5.25) + 
  geom_pointrange(alpha = 0.7) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) +
  labs(x = expression("DIP ("*mu*"M)"),
       y = expression("Symbiont Density ("~10^6*" cells/"~cm^-2*")"),
       color = element_blank()) + 
  scale_color_brewer(palette="Set2")
SymbP
grid.arrange(SymbN,SymbP,ncol=2)


symb3 <- symb2%>%mutate(NadjQuantBins = cut(Nadj, 
                                          breaks = unique(quantile(Nadj,probs=seq.int(0,1, by=1/4))), 
                                          include.lowest=TRUE))
symb3 <- symb3%>%mutate(PadjQuantBins = cut(Padj, 
                                          breaks = unique(quantile(Padj,probs=seq.int(0,1, by=1/4))), 
                                          include.lowest=TRUE))
N <- ggplot(symb3, aes(x = NadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Proportion of Reported Change in Symbiont Density",
       fill = "Symbiont Density") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(size = 15)) +
  scale_x_discrete(labels = c("[0.08,3]" = "0.08-3.0",
                              "(3,10]" = "3.0-10",
                              "(10,20]" = "10-20",
                              "(20,128]" = "20-128"))
P <- ggplot(symb3, aes(x = PadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "",
       fill = "Symbiont Density") +
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
  scale_x_discrete(labels = c("[0.02,0.05]" = "0.02-0.05",
                              "(0.05,0.228]" = "0.05-0.23",
                              "(0.228,2]" = "0.23-2.0"))
grid.arrange(N,P,ncol=2)


### dotplot
dotplotoverview <- ggplot(symb2, 
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
       size = "Treatment\nEffect Size", color = "Symbiont Density", 
       shape = "Optimal Density") + 
  guides(shape = guide_legend(order = 3, override.aes = list(size=4)),
         color = guide_legend(order = 2, override.aes = list(size=4)),
         size = guide_legend(order = 1)) +
  scale_shape_manual(values=c(17, 16))


### COVARIANCE MATRICES for meta-analysis
newlist <- list(NA)
for (i in seq(1,length(covar))) {
  newlist[i] <- list(covar[[i]]$S)
}

### LOOKING AT SMD
## N + P with zero
mod5 <- mixmeta(smd ~ 0 + Nadj + Padj, 
                random = ~1 | Reference, 
                data = symb2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod5) # I2 = 76.8

## N + P + N*P with zero
mod6 <- mixmeta(smd ~ 0 + Nadj + Padj + Nadj:Padj, 
                random = ~ 1| Reference, 
                data = symb2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod6) # I2 = 77

## log10N + log10P with zero
mod7 <- mixmeta(smd ~ 0 + log10N + log10P, 
                random = ~ 1 | Reference, 
                data = symb2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod7) # I2 = 73.2

## log10N + log10P + log10N*log10P with zero
mod8 <- mixmeta(smd ~ 0 + log10N + log10P + log10N:log10P, 
                random = ~ 1 | Reference, 
                data = symb2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod8) # I2 = 73.3

## N + P with simple random
mod9 <- mixmeta(smd ~ Nadj + Padj, 
                random = ~ 1 | Reference, 
                data = symb2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod9) # I2 = 74.4

## N + P + N*P with simple random
mod10 <- mixmeta(smd ~ Nadj + Padj + Nadj:Padj, 
                 random = ~ 1| Reference, 
                 data = symb2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod10) # I2 = 74.6

summary(symb2$NO3)
summary(symb2$NH4)
symb2$NO3adj <- ifelse(symb2$NO3 == 0, 0.1, symb2$NO3)
symb2$NH4adj <- ifelse(symb2$NH4 == 0, 0.1, symb2$NH4)

## log10N + log10P with simple random
mod11 <- mixmeta(smd ~ log10N + log10P, 
                 random = ~ 1 | Reference, 
                 data = symb2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod11) # I2 = 69.1%

mod11b <- mixmeta(smd ~ log10(NO3adj) + log10(NH4adj) + log10P, 
                 random = ~ 1 | Reference, 
                 data = symb2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod11b)
AIC(mod11, mod11b) # not an improvement

## log10N + log10 P + log10N*log10P with simple random
mod12 <- mixmeta(smd ~ log10N + log10P + log10N:log10P, 
                 random = ~ 1 | Reference, 
                 data = symb2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod12) # I2 = 69.9
AIC(mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12) # mod 11 & 12
BIC(mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12) # mod 11 & 12

## log10N + log10P with simple random
mod11b <- mixmeta(smd ~ log10N + log10P + Exposure.in.days, 
                 random = ~ 1 | Reference, 
                 data = symb2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod11b) 
AIC(mod11, mod11b)
BIC(mod11, mod11b)

resid <- resid(mod11)
fitted <- fitted(mod11)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)
## pretty good residuals 


### nonlinear models
## generating 1st & 3rd Q for nonlinear models
summary(symb2$log10N)
summary(symb2$log10P)
summary(symb2$Nadj)
summary(symb2$Padj)

## simple random
mod14 <- mixmeta(smd ~ 0 + ns(log10N, knots=c(0.4771,1.3010)) + 
                   ns(log10P, knots=c(-1.6990,-0.6437)), 
                 data=symb2,
                 random= ~ 1| Reference, 
                 method="ml", 
                 control=list(addSlist=newlist))
summary(mod14) # I2 = 70.6

mod14b <- mixmeta(smd ~ 0 + ns(Nadj, knots = c(3.00, 20.00)) +
                    ns(Padj, knots = c(0.0200, 0.2275)), 
                  data=symb2,
                  random= ~ 1| Reference, 
                  method="ml", 
                  control=list(addSlist=newlist))
summary(mod14b)

## with quartiles
mod14quart <- mixmeta(smd ~ 0+ ns(log10N, knots=c(0.4771,1.0000,1.3010)) + 
                        ns(log10P, knots=c(-1.6990,-1.3010,-0.6437)), 
                      data=symb2,
                      random= ~ 1|Reference, 
                      method="ml", 
                      control=list(addSlist=newlist))
summary(mod14quart) # I2 = 71.3

mod14quartb <- mixmeta(smd ~ 0+ ns(Nadj, knots=c(3.00,10.00,20)) + 
                         ns(Padj, knots=c(.02,.05,0.2275)), 
                       data=symb2,
                       random= ~ 1|Reference, 
                       method="ml", 
                       control=list(addSlist=newlist))
summary(mod14quartb)

AIC(mod11, mod14, mod14b, mod14quart, mod14quartb)
BIC(mod11, mod14, mod14b, mod14quart, mod14quartb)
## mod14b best with AIC, but mod11 and mod14b are equiv with BIC
summary(mod14b)
resid <- resid(mod14b)
fitted <- fitted(mod14b)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)
## VS. 
resid <- resid(mod11)
fitted <- fitted(mod11)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)
## residuals look pretty equivalent

### checking for benefit of adding exposure duration
mod14c <- mixmeta(smd ~ 0 +ns(Nadj, knots = c(3.00, 20.00)) +
                    ns(Padj, knots = c(0.0200, 0.2275)) + Exposure.in.days, 
                  data=symb2,
                  random= ~ 1| Reference, 
                  method="ml", 
                  control=list(addSlist=newlist))
summary(mod14c)

mod14c <- mixmeta(smd ~ 0 +ns(Nadj, knots = c(3.00, 20.00)) +
                    ns(Padj, knots = c(0.0200, 0.2275)) + Exposure.in.days, 
                  data=symb2,
                  random= ~ 1| Reference, 
                  method="ml", 
                  control=list(addSlist=newlist))
summary(mod14c)
AIC(mod14b, mod14c)
## mod14b is better
mod14d <- mixmeta(smd ~ 0 +ns(NO3adj, knots = c(0.100, 4.910)) + ns(NH4adj, knots = c(0.100, 15.000)) + 
                    ns(Padj, knots = c(0.0200, 0.2275)) + Exposure.in.days, 
                  data=symb2,
                  random= ~ 1| Reference, 
                  method="ml", 
                  control=list(addSlist=newlist))
summary(mod14d)


# PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS
## with linear model 3
newdata.predictN <- data.frame(log10N = symb2$log10N,
                               log10P = median(symb2$log10P))
newdata.predictP <- data.frame(log10N = median(symb2$log10N),
                               log10P = symb2$log10P)
predN <- predict(mod11, newdata=newdata.predictN, ci=TRUE)
predP <- predict(mod11, newdata=newdata.predictP, ci=TRUE)
testCIN <- cbind(symb2, predN)
testCIP <- cbind (symb2, predP)
min_conc_N <- testCIN %>% 
  mutate(overlap0 = 0 >= ci.lb) %>% #& 0 <= ci.ub) %>% 
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
        legend.position = c(0.15, 0.75),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) +  
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b")  +
  labs(color = expression("DIP ("*mu*"M)")) +
  geom_ribbon(aes(x = Nadj, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  geom_line(aes(x = Nadj, y = fit), inherit.aes=FALSE) +
  geom_vline(xintercept=3.4, linetype="dashed", color = "maroon") 
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
grid.arrange(ESNlin,ESPlin,ncol=2)
summary(mod11)


# PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS
## with nonlinear model 14b
newdata.predictN2 <- data.frame(Nadj = symb2$Nadj,
                               Padj = median(symb2$Padj))
newdata.predictP2 <- data.frame(Nadj = median(symb2$Nadj),
                                Padj = symb2$Padj)
predNnonlin <- predict(mod14b, newdata=newdata.predictN2, ci=TRUE)
predPnonlin <- predict(mod14b, newdata=newdata.predictP2, ci=TRUE)
testCINnonlin <- cbind(symb2, predNnonlin)
testCIPnonlin <- cbind (symb2, predPnonlin)
min_conc_Nnonlin <- testCINnonlin %>% 
  mutate(overlap0 = 0 >= ci.lb & 0 <= ci.ub) %>% 
  filter(overlap0==FALSE) %>% 
  summarize(min_N = min(Nadj),
            max_N = max(Nadj),
            min_P = min(Padj),
            max_P = max(Padj))

## and now plotting for nitrogen as a predictor
ESNnonlin <- ggplot(testCINnonlin, 
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
        legend.position = c(0.15, 0.75),
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
ESNnonlin

ESPnonlin <- ggplot(testCIPnonlin, 
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
ESPnonlin
grid.arrange(ESNnonlin,ESPnonlin,ncol=2)


