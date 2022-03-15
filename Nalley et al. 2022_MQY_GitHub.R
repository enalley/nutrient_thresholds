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
allmqy <- as.data.frame(read.csv("ALL-MQY.csv", header=TRUE))
str(allmqy) ## everything is factor vs. numeric

# creating new variable that has both Genus and species
allmqy <- allmqy %>%
  mutate(Gsp = paste(Genus, Species, sep = "_")) 
head(allmqy)

### estimating ambient concentrations of N & P if listed as zero
allmqy$Nadj <- ifelse(allmqy$TOTAL.N == 0, 0.1, allmqy$TOTAL.N)
allmqy$Padj <- ifelse(allmqy$PO4 == 0, 0.02, allmqy$PO4)
### adding log variables
allmqy$log10N <- log10(allmqy$Nadj)
allmqy$log10P <- log10(allmqy$Padj)

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

MQYN <- ggplot(allmqy, 
               aes(x = Nadj, y = Response.level, 
                   ymin = Response.level-SE,
                   ymax = Response.level+SE,
                   color = Control)) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +  ylim(0,1) +
  geom_smooth(method = "loess", color = "darkgrey", fill = "lightgrey") +
  geom_pointrange(alpha = 0.7) + theme_bw() +
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
       y = "Max. Photo. Efficiency ("~F[v]*"/"~F[m]*")",
       color =  element_blank()) +
  scale_color_brewer(palette="Set2")
MQYN

MQYP <- ggplot(allmqy, 
               aes(x = Padj, y = Response.level, 
                   ymin = Response.level-SE,
                   ymax = Response.level+SE,
                   color = Control)) + 
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") + ylim(0,1) +
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
       y = "Max. Photo. Efficiency ("~F[v]*"/"~F[m]*")",
       color =  element_blank()) +
  scale_color_brewer(palette="Set2")
MQYP
grid.arrange(MQYN,MQYP,ncol=2)


allmqy3 <- allmqy2%>%mutate(NadjQuantBins = cut(Nadj, 
                                                breaks = unique(quantile(Nadj,probs=seq.int(0,1, by=1/4))), 
                                                include.lowest=TRUE))
allmqy3 <- allmqy3%>%mutate(PadjQuantBins = cut(Padj, 
                                                breaks = unique(quantile(Padj,probs=seq.int(0,1, by=1/4))), 
                                                include.lowest=TRUE))
N <- ggplot(allmqy3, aes(x = NadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Proportion of Reported Change in Photo. Efficiency",
       fill = "") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(size = 15)) +
  scale_x_discrete(labels = c("[0.33,2.3]" = "0.3-2.3",
                              "(2.3,5.15]" = "2.3-5.2",
                              "(5.15,8.95]" = "5.2-9.0",
                              "(8.95,128]" = "9.0-128"))
P <- ggplot(allmqy3, aes(x = PadjQuantBins, fill = smdpos)) +
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
  scale_x_discrete(labels = c("[0.02,0.0333]" = "0.02-0.03",
                              "(0.0333,0.294]" = "0.03-0.3",
                              "(0.294,0.645]" = "0.3-0.6",
                              "(0.645,0.7]" = "0.6-0.7"))
grid.arrange(N,P,ncol=2)



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
## N + P with zero
mod5 <- mixmeta(smd ~ 0 + Nadj + Padj, 
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod5) # I2 = 73.5

## N + P + N*P with zero
mod6 <- mixmeta(smd ~ 0 + Nadj + Padj + Nadj:Padj, 
                random = ~ 1| Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod6) 
mod6b <- mixmeta(smd ~ 0 + Nadj + Padj + Nadj:Padj + Exposure.in.days, 
                random = ~ 1| Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod6b)
AIC(mod6, mod6b) #mod6 is better
BIC(mod6, mod6b)


## log10N + log10P with zero
mod7 <- mixmeta(smd ~ 0 + log10N + log10P, 
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod7) # I2 = 77.5

## log10N + log10P + log10N*log10P with zero
mod8 <- mixmeta(smd ~ 0 + log10N + log10P + log10N:log10P, 
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod8) # I2 = 78.9

## N + P with simple random
mod9 <- mixmeta(smd ~ Nadj + Padj, 
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod9) # I2 = 72.5

mod9b <- mixmeta(smd ~ Nadj + Padj + Exposure.in.days, 
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod9b) 

mod9c <- mixmeta(smd ~ Nadj + Padj + Study, 
                random = ~ 1 | Reference, 
                data = allmqy2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod9c) 
AIC(mod9, mod9b, mod9c)

## N + P + N*P with simple random
mod10 <- mixmeta(smd ~ Nadj + Padj + Nadj:Padj, 
                 random = ~ 1| Reference, 
                 data = allmqy2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod10) # I2 = 70.4

## log10N + log10P with simple random
mod11 <- mixmeta(smd ~ log10N + log10P, 
                 random = ~ 1 | Reference, 
                 data = allmqy2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod11) # I2 = 77.5

## log10N + log10 P + log10N*log10P with simple random
mod12 <- mixmeta(smd ~ log10N + log10P + log10N:log10P, 
                 random = ~ 1 | Reference, 
                 data = allmqy2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod12) # I2 = 72.1
AIC(mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12) # mod6
BIC(mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod12) # mod6
summary(mod6)
resid <- resid(mod6)
fitted <- fitted(mod6)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)
## pretty good residuals

## generating 1st & 3rd Q for nonlinear models
summary(allmqy2$log10N)
summary(allmqy2$log10P)
summary(allmqy2$Nadj)
summary(allmqy2$Padj)

## nonlin
## simple random
mod14 <- mixmeta(smd ~ 0 + ns(log10N, knots = c(0.3617, 0.9409)) +
                   ns(log10P, knots = c(-1.5584, -0.1942)),
                 data=allmqy2,
                 random= ~ 1| Reference, 
                 method="ml", 
                 control=list(addSlist=newlist))
summary(mod14) 

mod14b <- mixmeta(smd ~ 0 +  ns(Nadj, knots = c(2.30, 8.95)) +
                    ns(Padj, knots = c(0.03325, 0.64550)),
                  data=allmqy2,
                  random= ~ 1| Reference, 
                  method="ml", 
                  control=list(addSlist=newlist))
summary(mod14b)

summary(allmqy2$log10N)
summary(allmqy2$log10P)
## with quartiles
mod14quart <- mixmeta(smd ~ 0 + ns(log10N, knots=c(0.3617,0.7083,0.9409)) + 
                        ns(log10P, knots=c(-1.5584,-0.5317,-0.1942)) , 
                      data=allmqy2,
                      random= ~ 1 |Reference, 
                      method="ml", 
                      control=list(addSlist=newlist))
summary(mod14quart) 

mod14quartb <- mixmeta(smd ~ 0+  ns(Nadj, knots=c(2.30,5.15,8.95)) + 
                         ns(Padj, knots=c(0.03325,0.29400,0.64550)), 
                       data=allmqy2,
                       random= ~ 1|Reference, 
                       method="ml", 
                       control=list(addSlist=newlist))
summary(mod14quartb)  

AIC(mod6, mod9, mod14, mod14b, mod14quart, mod14quartb)
BIC(mod6, mod9, mod14, mod14b, mod14quart, mod14quartb)
## mod 6 is still the best
summary(mod6)

resid <- resid(mod14quart)
fitted <- fitted(mod14quart)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)

resid <- resid(mod9)
fitted <- fitted(mod9)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)
AIC(mod9, mod14quart, mod14)
BIC(mod9, mod14quart, mod14)
## mod14qurt is better, but mod14 is not sig worse in AIC & is simpler

# PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS
## LINEAR - MOD9 - best fit model without interaction term
summary(mod9)
newdata.predictN <- data.frame(Nadj = allmqy2$Nadj,
                               Padj = median(allmqy2$Padj))
newdata.predictP <- data.frame(Nadj = median(allmqy2$Nadj),
                               Padj = allmqy2$Padj)
predN <- predict(mod9, newdata=newdata.predictN, ci=TRUE)
predP <- predict(mod9, newdata=newdata.predictP, ci=TRUE)
testCIN <- cbind(allmqy2, predN)
testCIP <- cbind(allmqy2, predP)

## and now plotting for nitrogen as a predictor
EffN <- ggplot(testCIN, 
               aes(x = Nadj,
                   y = smd, color = Padj,
                   ymin = smd-vmd,
                   ymax = smd+vmd)) + 
  geom_pointrange(size = 1, alpha = 0.7) +
#  ylim(-5,5) +
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

