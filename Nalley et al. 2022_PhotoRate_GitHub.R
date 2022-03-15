#### NALLEY ET AL. 2022 NUTRIENT THRESHOLD ANALYSIS #####
### Photosynthetic Rate (µmol O2 cm-2 day-1) ###

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

allpho <- as.data.frame(read.csv("ALL-PhotoSA.csv", header=TRUE))
str(allpho) ## everything is factor vs. numeric
allpho <- allpho[,-c(19:27)]

# creating new variable that has both Genus and species
allpho <- allpho %>%
  mutate(Gsp = paste(Genus, Species, sep = "_")) 
head(allpho)

### estimating ambient concentrations of N & P if listed as zero
allpho$Nadj <- ifelse(allpho$TOTAL.N == 0, 0.1, allpho$TOTAL.N)
allpho$Padj <- ifelse(allpho$PO4 == 0, 0.02, allpho$PO4)
### adding log variables
allpho$log10N <- log10(allpho$Nadj)
allpho$log10P <- log10(allpho$Padj)

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

PhotoN <- ggplot(allpho, 
                 aes(x = Nadj, y = Response.level, 
                     ymin = Response.level-SE,
                     ymax = Response.level+SE,
                     color = Control)) +
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") + ylim(-1,9) + 
  geom_smooth(method = "loess", color = "darkgrey", fill = "lightgrey") +
  geom_pointrange(alpha = 0.7) + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = c(0.2, 0.85),
        legend.background = element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=12)) +
  labs(x =expression("DIN ("*mu*"M)"),
       y = expression("Photo. Rate ("*mu*"mol"~O[2]*""~cm^-2*""~day^-1*")"),
       color = element_blank()) +
  scale_color_brewer(palette="Set2")
PhotoN

PhotoP <- ggplot(allpho, 
                 aes(x = Padj, y = Response.level, 
                     ymin = Response.level-SE,
                     ymax = Response.level+SE,
                     color = Control)) +
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +  
  geom_smooth(method = "loess", color = "darkgrey", fill = "lightgrey") +
  geom_pointrange(alpha = 0.7) + ylim(-1,9) + theme_bw() +
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
       y = expression("Photo. Rate ("*mu*"mol"~O[2]*""~cm^-2*""~day^-1*")"),
       color = element_blank()) +
  scale_color_brewer(palette="Set2")
grid.arrange(PhotoN,PhotoP,ncol=2)



allpho3 <- allpho2%>%mutate(NadjQuantBins = cut(Nadj, 
                                            breaks = unique(quantile(Nadj,probs=seq.int(0,1, by=1/4))), 
                                            include.lowest=TRUE))
allpho3 <- allpho3%>%mutate(PadjQuantBins = cut(Padj, 
                                            breaks = unique(quantile(Padj,probs=seq.int(0,1, by=1/4))), 
                                            include.lowest=TRUE))

N <- ggplot(allpho3, aes(x = NadjQuantBins, fill = smdpos)) +
  geom_bar(position = "fill") +
  labs(x = expression("DIN ("*mu*"M)"),
       y = "Proportion of Reported Change in Photo. Rate",
       fill = "") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(size = 15)) +
  scale_x_discrete(labels = c("[0.1,0.2]" = "0.1-0.2",
                              "(0.2,2.97]" = "0.2-3.0",
                              "(2.97,5]" = "3.0-5.0",
                              "(5,39.1]" = "5.0-39.1"))
P <- ggplot(allpho3, aes(x = PadjQuantBins, fill = smdpos)) +
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
  scale_x_discrete(labels = c("[0.02,0.05]" = "0.02-0.05",
                              "(0.05,0.2]" = "0.05-0.2",
                              "(0.2,0.7]" = "0.2-0.7",
                              "(0.7,5.14]" = "0.7-5.14"))
grid.arrange(N,P,ncol=2)

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
## N + P with zero
mod5 <- mixmeta(smd ~ 0 + Nadj + Padj, 
                random = ~ 1 | Reference, 
                data = allpho2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod5) 

## N + P + N*P with zero
mod6 <- mixmeta(smd ~ 0 + Nadj + Padj + Nadj:Padj, 
                random = ~ 1| Reference, 
                data = allpho2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod6) 

## log10N + log10P with zero
mod7 <- mixmeta(smd ~ 0 + log10N + log10P, 
                random = ~ 1 | Reference, 
                data = allpho2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod7) 

## log10N + log10P + log10N*log10P with zero
mod8 <- mixmeta(smd ~ 0 + log10N + log10P + log10N*log10P, 
                random = ~ 1 | Reference, 
                data = allpho2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod8) 

## N + P with simple random
mod9 <- mixmeta(smd ~ Nadj + Padj, 
                random = ~ 1 | Reference, 
                data = allpho2, method = "ml", 
                control = list(addSlist = newlist))
summary(mod9) 

## N + P + N*P with simple random
mod10 <- mixmeta(smd ~ Nadj + Padj + Nadj:Padj, 
                 random = ~ 1| Reference, 
                 data = allpho2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod10) 

## log10N + log10P with simple random
mod11 <- mixmeta(smd ~ log10N + log10P, 
                 random = ~ 1 | Reference, 
                 data = allpho2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod11) 

mod11b <- mixmeta(smd ~ 0 + log10N + log10P, 
                 random = ~ 1 | Reference, 
                 data = allpho2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod11b)
AIC(mod11, mod11b)

## log10N + log10 P + log10N*log10P with simple random
mod12 <- mixmeta(smd ~ log10N + log10P + log10N:log10P, 
                 random = ~ 1 | Reference, 
                 data = allpho2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod12) # I2 = 48.4
AIC(mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod11b, mod12) 
BIC(mod5, mod6, mod7, mod8, mod9, mod10, mod11, mod11b, mod12) 
summary(mod8)
resid <- resid(mod8)
fitted <- fitted(mod8)
plot(fitted, resid)
abline(0,0)
qqnorm(resid)
qqline(resid)
## pretty good residuals 

## generating 1st & 3rd Q for nonlinear models
summary(allpho2$log10N)
summary(allpho2$log10P)
summary(allpho2$Nadj)
summary(allpho2$Padj)

## nonlin
mod13 <- mixmeta(smd ~ 0 + ns(log10N, knots = c(-0.6990, 0.6990)) +
                   ns(log10P, knots = c(-1.3010, -0.1664)),
                 random = ~ 1 | Reference, 
                 data = allpho2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod13) # I2 = 59.2

mod13b <- mixmeta(smd ~ 0 + ns(Nadj, knots = c(0.200, 5.000)) +
                    ns(Padj, knots = c(0.050, 0.700)),
                  random = ~ 1 | Reference, 
                  data = allpho2, method = "ml", 
                  control = list(addSlist = newlist))
summary(mod13b) # I2 = 59.2
AIC(mod11, mod13, mod13b)

mod13c <- mixmeta(smd ~ 0 + ns(Nadj, knots = c(0.200, 5.000)) +
                    ns(Padj, knots = c(0.050, 0.700)) + Exposure.in.days, 
                  random = ~ 1 | Reference, 
                  data = allpho2, method = "ml", 
                  control = list(addSlist = newlist))
summary(mod13c)

mod13d <- mixmeta(smd ~ 0 + ns(log10N, knots = c(-0.6990, 0.6990)) +
                   ns(log10P, knots = c(-1.3010, -0.1664)) + Exposure.in.days, 
                 random = ~ 1 | Reference, 
                 data = allpho2, method = "ml", 
                 control = list(addSlist = newlist))
summary(mod13d) 

## with quartiles
mod13quart <- mixmeta(smd ~ 0 + ns(Nadj, knots=c(0.200,2.970,5.000)) + 
                        ns(Padj, knots=c(0.050, 0.200, 0.70)), 
                      data=allpho2,
                      random= ~ 1|Reference, 
                      method="ml", 
                      control=list(addSlist=newlist))
summary(mod13quart)

## with quartiles
mod14quart <- mixmeta(smd ~ 0 + ns(log10N, knots=c(-0.6990,0.4727,0.6990)) + 
                        ns(log10P, knots=c(-1.3010,-0.6990,-0.1664)), 
                      data=allpho2,
                      random= ~ 1 | Reference, 
                      method="ml", 
                      control=list(addSlist=newlist))
summary(mod14quart)

AIC(mod8, mod13, mod13b, mod13c, mod13d, mod13quart, mod14quart)
BIC(mod8, mod13, mod13b, mod13c, mod13d, mod13quart, mod14quart)
## mod8 is best
summary(mod8)

### plotting with predictions
newdata.predictN <- data.frame(log10N = allpho2$log10N,
                               log10P = median(allpho2$log10P))
newdata.predictP <- data.frame(log10N = median(allpho2$log10N),
                               log10P = allpho2$log10P)
predN <- predict(mod8, newdata=newdata.predictN, ci=TRUE, vcov = TRUE)
predP <- predict(mod8, newdata=newdata.predictP, ci=TRUE, vcov = TRUE)
testCIN <- cbind(allpho2, predN)
testCIP <- cbind (allpho2, predP)

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
        legend.position = c(0.15, 0.2),
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
        legend.position = c(0.85, 0.2),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.background = element_blank()) +  
  scale_x_log10(breaks = breaks_log10) +  
  annotation_logticks(sides = "b") +
  labs(color = expression("DIN ("*mu*"M)"))  
grid.arrange(ESNlin,ESPlin,ncol=2)


