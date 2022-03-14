#### NALLEY ET AL. 2022 NUTRIENT STRESSOR THRESHOLDS ####
#### ChlA Concentration (ug chl/cm^2) ####
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
chl <- as.data.frame(read.csv("ALL-ChlA_SA.csv", header=TRUE))
str(chl) ## everything is factor vs. numeric
chl <- chl[,c(1:17)]
str(chl)

# creating new variable that has both Genus and species
chl <- chl %>%
  mutate(Gsp = paste(Genus, Species, sep = "_")) 
head(chl)
### estimating ambient concentrations of N & P if listed as zero
chl$Nadj <- ifelse(chl$TOTAL.N == 0, 0.1, chl$TOTAL.N)
chl$Padj <- ifelse(chl$PO4 == 0, 0.02, chl$PO4)
### adding log variables
chl$log10N <- log10(chl$Nadj)
chl$log10P <- log10(chl$Padj)

### ordering data before calculating hedges
str(chl)
chl.ord <- chl[,c(1,3:ncol(chl))]
chl.ord <- chl[order(chl$Reference, chl$Control),]
### calculating hedges G (std. mean difference)
### Reference = experiment (within studies)
covar <- by(chl.ord, chl.ord$Reference, function(x) 
  covar.smd(Response.level, SD, N, "smd", method="hedges", data = x))
chl.ord$smd <- unlist(lapply(covar, function(x) x$y))
chl.ord$vmd <- unlist(lapply(covar, function(x) x$v))
str(chl.ord) 
test <- chl.ord[,c("Reference", "Control", "smd")]

### removing controls for modeling
chl2 <- subset(chl.ord, Control=="exp")
str(chl2)
### adding variable for inc/dec chla
chl2$smdpos <- as.factor(ifelse(chl2$smd > 0, "Increase", "Decrease"))
str(chl2)

## Summary of data 
str(chl2)
## 22 papers with 37 experiments within
ChlTreatSummary <- chl2 %>% 
  summarize(Nmin = min(Nadj),
            Nmax = max(Nadj),
            Pmin = min(Padj),
            Pmax = max(Padj),
            DurMean = mean(Exposure.in.days),
            DurMin = min(Exposure.in.days),
            DurMax = max(Exposure.in.days))


ChlN <- ggplot(chl, 
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
        legend.position = c(0.15, 0.8),
        legend.text=element_text(size=16),
        legend.title=element_blank()) + 
  labs(x = expression("DIN ("*mu*"M)"),
       y = expression("Chl-a ("*mu*"g"~cm^-2*")"),
       color = element_blank()) +
  scale_color_brewer(palette="Set2")
ChlN

ChlP <- ggplot(chl, 
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
        text = element_text(size=20),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) + 
  labs(x = expression("DIP ("*mu*"M)"),
       y = expression("Chl-a ("*mu*"g"~cm^-2*")"),
       color = expression("Total N "*mu*"M")) +
  scale_color_brewer(palette="Set2")
ChlP
grid.arrange(ChlN,ChlP,ncol=2)

## dot plot
ggplot(chl2, 
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
       size = "Treatment\nEffect Size", color = "Chl-a Conc.") + 
  guides(color = guide_legend(order = 2, override.aes = list(size=4)),
         size = guide_legend(order = 1))

chl3 <- chl2%>%mutate(NadjQuantBins = cut(Nadj, 
                                             breaks = unique(quantile(Nadj,probs=seq.int(0,1, by=1/4))), 
                                             include.lowest=TRUE))
chl3 <- chl3%>%mutate(PadjQuantBins = cut(Padj, 
                                          breaks = unique(quantile(Padj,probs=seq.int(0,1, by=1/4))), 
                                          include.lowest=TRUE))
summary(chl3$Padj)
N <- ggplot(chl3, aes(x = NadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Proportion of Reported Change in Chl-a Conc.",
       fill = "Chl-a Conc.") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(size = 15)) +
  scale_x_discrete(labels = c("[0.1,1.9]" = "0.1-1.9",
                              "(1.9,5.09]" = "1.9-5.1",
                              "(5.09,20]" = "5.1-20",
                              "(20,50]" = "20-50"))
P <- ggplot(chl3, aes(x = PadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIP ("*mu*"M)"),
       y = "Proportion",
       fill = "Chl-a Conc.") +
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
  scale_x_discrete(labels = c("[0.02,0.12]" = "0.02-0.1",
                              "(0.12,0.3]" = "0.1-0.3",
                              "(0.3,5.14]" = "0.3-5.1"))
grid.arrange(N,P,ncol=2)



### COVARIANCE MATRICES for meta-analysis
newlist <- list(NA)
for (i in seq(1,length(covar))) {
  newlist[i] <- list(covar[[i]]$S)
}

### LOOKING AT SMD
## N + P with zero
mod5 <- mixmeta(smd ~ 0+ Nadj + Padj, 
                random = ~ 1 | Reference, 
                data = chl2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod5) # I2 = 52

## N + P + N*P with zero
mod6 <- mixmeta(smd ~ 0+Nadj + Padj + Nadj:Padj, 
                random = ~ 1| Reference, 
                data = chl2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod6) # I2 = 50

## log10N + log10P with zero
mod7 <- mixmeta(smd ~ 0+ log10N + log10P, 
                random = ~ 1 | Reference, 
                data = chl2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod7) # I2 = 46

## log10N + log10P + log10N*log10P with zero
mod8 <- mixmeta(smd ~ 0+log10N + log10P + log10N:log10P, 
                random = ~ 1 | Reference, 
                data = chl2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod8) # I2 = 47

## N + P with simple random
mod9 <- mixmeta(smd ~ Nadj + Padj, 
                random = ~ 1 | Reference, 
                data = chl2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod9) # I2 = 46

## N + P + N*P with simple random
mod10 <- mixmeta(smd ~ Nadj + Padj + Nadj:Padj, 
                 random = ~ 1| Reference, 
                 data = chl2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod10) # I2 = 47

## log10N + log10P with simple random
mod11 <- mixmeta(smd ~ log10N + log10P, 
                 random = ~ 1 | Reference, 
                 data = chl2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod11) # I2 = 43

mod11b <- mixmeta(smd ~ log10N + log10P, 
                 random = ~ 1 | Reference, 
                 data = chl2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod11b)

## log10N + log10 P + log10N*log10P with simple random
mod12 <- mixmeta(smd ~ log10N + log10P + log10N:log10P, 
                 random = ~ 1 | Reference, 
                 data = chl2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod12) # I2 = 41
AIC(mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12) # mod 11 & 12
BIC(mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12) # mod 11 & 12

## generating 1st & 3rd Q for nonlinear models
summary(chl2$Nadj)
summary(chl2$Padj)
summary(chl2$log10N)
summary(chl2$log10P)

## simple random
mod14 <- mixmeta(smd ~ 0 + ns(log10N, knots=c(0.2788,1.3010)) + 
                   ns(log10P, knots=c(-1.6990,-0.5229)), 
                 data=chl2,
                 random= ~ 1| Reference, 
                 method="ml", 
                 control=list(addSlist=newlist))
summary(mod14) # I2 = 41

mod14b <- mixmeta(smd ~ 0 + ns(Nadj, knots = c(1.90, 20.00)) +
                    ns(Padj, knots = c(0.0200, 0.3000)), 
                  data=chl2,
                  random= ~ 1| Reference, 
                  method="ml", 
                  control=list(addSlist=newlist))
summary(mod14b) ## I2 = 30

mod14c <- mixmeta(smd ~ 0 + ns(Nadj, knots = c(1.90, 20.00)) +
                    ns(Padj, knots = c(0.0200, 0.3000)) + Exposure.in.days, 
                  data=chl2,
                  random= ~ 1| Reference, 
                  method="ml", 
                  control=list(addSlist=newlist))
summary(mod14c)
AIC(mod14b, mod14c)

## with quartiles
mod14quart <- mixmeta(smd ~ 0+ ns(log10N, knots=c(0.2788,0.7067,1.3010)) + 
                        ns(log10P, knots=c(-1.6990,-0.9208,-0.5229)), 
                      data=chl2,
                      random= ~ 1|Reference, 
                      method="ml", 
                      control=list(addSlist=newlist))
summary(mod14quart) # I2 = 27

mod14quartb <- mixmeta(smd ~ 0+ ns(Nadj, knots=c(1.90,5.09,20)) + 
                         ns(Padj, knots=c(.02,0.1200,0.3000)), 
                       data=chl2,
                       random= ~ 1|Reference, 
                       method="ml", 
                       control=list(addSlist=newlist))
summary(mod14quartb) # I2 = 29
AIC(mod11, mod14, mod14b, mod14quart, mod14quartb)
BIC(mod11, mod14, mod14b, mod14quart, mod14quartb)
## mod14b best with both

summary(mod14b)
resid <- resid(mod14b)
fitted <- fitted(mod14b)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)


# PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS
newdata.predictN <- data.frame(Nadj = chl2$Nadj,
                               Padj = median(chl2$Padj))
newdata.predictP <- data.frame(Nadj = median(chl2$Nadj),
                               Padj = chl2$Padj)
predN <- predict(mod14b, newdata=newdata.predictN, ci=TRUE)
predP <- predict(mod14b, newdata=newdata.predictP, ci=TRUE)
testCIN <- cbind(chl2, predN)
testCIP <- cbind(chl2, predP)
min_conc_d <- testCIN %>% 
  mutate(overlap0 = 0 >= ci.lb & 0 <= ci.ub) %>% 
  filter(overlap0==FALSE) %>% 
  summarize(min_N = min(Nadj),
            max_N = max(Nadj),
            min_P = min(Padj),
            max_P = max(Padj))

## and now plotting for nitrogen as a predictor
EffN <- ggplot(testCIN, 
       aes(x =Nadj,
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
        legend.position = c(0.85, 0.85),
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

