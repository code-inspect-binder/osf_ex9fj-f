####### DeFiEn Project - Analysis to accompany the manuscript #######
# To submit with DeFiEn manuscript to SSLA

#clear and set up working directory
rm(list = ls())
# load required packages, functions, and lists (provided)
library(stopwords)
library(dplyr)
library(ggplot2)
library(psych)
library(lmerTest)
library(lme4)
library(tidyverse)
library(merTools)
stopwords <- stopwords() 
source("source.R")

load("data.rda")
#this archive file contains three dataframe: 
#joint (eye-movement to all interest areas); 
#joint.fix (eye-movements to fixated interest areas only)
#joint.fix.freq (eye-movements to fixated interest areas which have frequency values in the SUBTLEX corpus)
load("de_en_fi.rda") #joint.fix.freq broken down by language into three dataframes.
load("extremes.rda") #median and top values of three tests, with all combinations
load("ReadRateData.rda") #reading rate data with individual differences tests, in three dataframes broken down by language
load("joint.read_rate.rda") #a summary of reading rate (wmp) by participant

####### Descriptive statistics #######

### Descriptive statistics of the three component skills tests
### Table 1
parts = unique(joint[, c("Subject", "ART", "SpellTest", "Language", "Lextale")])
describeBy(parts[, c("ART", "SpellTest", "Lextale")], group = parts$Language, mat = T)[, c("group1", 'mean', 'median', "sd", "range")]
# 'NA' warning due to the fact that the English cohort did not complete the vocabulary test

### Multiple-regression mixed-effects models fitted separately to ART, spelling, and vocabulary scores with language as a predictor 
### Shows the significant differences between the cohorts as shown in Appendix A, Table A1, A2, and A3

### Table A1 - Author Recognition Test
art_lm = lm(ART ~ Language, data = parts)
summary(art_lm)

### Table A2 - Spelling Test
spell_lm = lm(SpellTest ~ Language, data = parts) 
summary(spell_lm)

### Table A3 - Vocabulary Test
lex_lm = lm(Lextale ~ Language, data = parts)
summary(lex_lm)

### boxplots of the individual differences tests
### Figure 1
par(mfrow = c(2,2))
boxplot(ART ~ Language, data = parts, main = "Author recognition test")
boxplot(SpellTest ~ Language, data = parts, main = "Spelling")
boxplot(Lextale ~ Language, data = parts, main = "Vocabulary")
par(mfrow = c(1,1)) 

### descriptive statistics for the eye-tracking variables of interest
### ToFD and others
# Table 2 part 1
describeBy(joint.fix$ToFD, group=joint.fix$Language, mat = T)[, c('group1', 'mean', 'sd', 'median', 'range')]
### Reading Rate
# Table 2 part 2
describe(fi.read_rate$ReadRate)[, c('mean', 'sd', 'median', 'range')] 
describe(de.read_rate$ReadRate)[, c('mean', 'sd', 'median', 'range')] 
describe(eng.read_rate$ReadRate)[, c('mean', 'sd', 'median', 'range')]

### descriptive statistics for the word characteristics of the sentence stimuli (length and frequency)
words = unique(joint.fix.freq[, c("length", "SUBTLEX", "Language")])
describeBy(words[, c("length", "SUBTLEX")], group = words$Language, mat = T)[, c('group1', 'mean', 'sd', 'median', 'range')]

####### Reading Rate #######
### Model that shows significant differences between the three languages in RR
# Reading rate as a function of language (German baseline)
# Table 3
mod.read_rate = lm(ReadRate ~ Language, data = joint.read_rate)
summary(mod.read_rate)

#boxplot visualizing the differences between the three languages
par(mfrow = c(1,1))
boxplot(ReadRate ~ Language, data = joint.read_rate, main = "Reading rate, words/minute")

####### Statistical Prediction for Reading Rate #######
### Baseline models for both German and Finnish
# Baseline model for German
# Reading rate as a function of component skills
# Table 4
read_rate.de = lmer(ReadRate ~ (SpellTest) + (ART) + (Lextale) + (1 |Subject) + (1 | text), data = de.read_rate)
summary(read_rate.de)

# Baseline model for Finnish
# Reading rate as a function of component skills
# Table 5
read_rate.fi = lmer(ReadRate ~ (SpellTest) + (ART) + (Lextale) +(1 |Subject) + (1 | text), data = fi.read_rate)
summary(read_rate.fi)

### Predicted reading rates for hypothetical German and Finish readers
# First, we used the predict function on the previous baseline models

## German
# Using the extremes.rda file (which has the top scores and median scores for all predictor variables)
# Top skills (last row of dataframe extremes):
predict(read_rate.de, newdata = extremes[8,], re.form = NA, allow.new.levels = T) -> pred.de
de.read_rate -> dat.de
dat.de$fitted = pred.de
# Median skills (first row of dataframe extremes):
predict(read_rate.de, newdata = extremes[1,], re.form = NA, allow.new.levels = T) -> pred.de.med
dat.de$fitted.med = pred.de.med

## Finnish
# Top skills:
predict(read_rate.fi, newdata = extremes[8,], re.form = NA, allow.new.levels = T) -> pred.fi
fi.read_rate -> dat.fi
dat.fi$fitted = pred.fi
# Median skills:
predict(read_rate.fi, newdata = extremes[1,], re.form = NA, allow.new.levels = T) -> pred.fi.med
dat.fi$fitted.med = pred.fi.med

## Predict functions for both the top and median skills
predFun <- function(fit) {
  predict(fit,extremes[8,], re.form = NA)
}
predFun.med <- function(fit) {
  predict(fit,extremes[1,], re.form = NA)
}

## Bootstrapping in 1000 iterations
# Note these values will differ slightly from the manuscript as it is a random estimation
# German Top:
bb.de <- bootMer(read_rate.de,nsim=1000,FUN=predFun, use.u = FALSE)
predMat.de <- bb.de$t
as.vector(quantile(predMat.de[,1], c(0.1,0.50,0.9))) -> q.de
# German Median:
bb.de.med <- bootMer(read_rate.de,nsim=1000,FUN=predFun.med, use.u = FALSE)
predMat.de.med <- bb.de.med$t
as.vector(quantile(predMat.de.med[,1], c(0.1,0.50,0.9))) -> q.de.med
# Finnish Top:
bb.fi <- bootMer(read_rate.fi,nsim=1000,FUN=predFun, use.u = FALSE)
predMat.fi <- bb.fi$t
as.vector(quantile(predMat.fi[,1], c(0.1,0.50,0.9))) -> q.fi
# Finnish Median:
bb.fi.med <- bootMer(read_rate.fi,nsim=1000,FUN=predFun.med, use.u = FALSE)
predMat.fi.med <- bb.fi.med$t
as.vector(quantile(predMat.fi.med[,1], c(0.1,0.50,0.9))) -> q.fi.med

### Values saved into file 'all'
data.frame(rbind(
  cbind(1, c(q.de, q.de.med)),
  cbind(2, c(as.vector(quantile(eng.read_rate$ReadRate, 0.1)), median(eng.read_rate$ReadRate), 
             as.vector(quantile(eng.read_rate$ReadRate, 0.9)))),
  cbind(3, c(q.fi, q.fi.med))
)) -> all

### Plot the values saved into 'all'
# Figure 2
par(mfrow = c(1,1))
plot(all[,1], all[, 2], xaxt = "n", pch = 19, xlab = 'Language', ylab = 'Reading Rate', type = "n") 
points(x = 1, y = q.de[2], pch = 18, col="grey60")
segments(x0 = 1, x1 = 1, y0 = q.de[1], y1 = q.de[3], col="grey60")
points(x = 1.1, y = q.de.med[2], pch = 18, col="grey60")
segments(x0 = 1.1, x1 = 1.1, y0 = q.de.med[1], y1 = q.de.med[3], col="grey60")
points(x = 3, y = q.fi[2], pch = 20, col="grey60")
segments(x0 = 3, x1 = 3, y0 = q.fi[1], y1 = q.fi[3], col="grey60")
points(x = 2.9, y = q.fi.med[2], pch = 20, col="grey60")
segments(x0 = 2.9, x1 = 2.9, y0 = q.fi.med[1], y1 = q.fi.med[3], col="grey60")
axis(side = 1, labels = c("German", "English", "Finnish"), at = 1:3, tick = T)
points(2.8, median(fi.read_rate$ReadRate), pch = 1)
points(1.2, median(de.read_rate$ReadRate), pch = 5)
points(2, 193, pch = 17)
points(2, 313, pch = 17)
points(2, 477, pch = 17)
legend("topright", legend = c("observed", "predicted"), pch = c(15, 15), col = c("black", "grey60"), bty = "n")
text(2, 193, "10%", pos = 4)
text(2, 313, "median", pos = 4)
text(2, 477, "90%", pos = 4)
text(1.2, median(de.read_rate$ReadRate), "median", pos = 4)
text(2.8, median(fi.read_rate$ReadRate), "median", pos = 2)
text(1, q.de[3], "top skills", pos = 4)
text(1.1, q.de.med[3], "median skills", pos = 4)
text(3, q.fi[3], "top skills", pos = 2)
text(2.9, q.fi.med[3], "median skills", pos = 2)
  
####### Total Fixation Duration #######

### Simple model of language and ToFD with by subject and by item intercepts
# Table A4 in Supplementary Materials
tofd_mod1 = lmer(log(ToFD) ~ Language + (1 | TargetWord) + (1 | Subject), data = joint.fix)
summary(tofd_mod1)

### Plot of ToFD and word length for each language -- split by skill level
# 'summarySE' is a function created in source.R
expr = summarySE(data = joint.fix, measurevar = "ToFD", groupvars = c("Language", "Skill", "length_bin"))
expr
# identify the levels of language
levels(expr$Language)[levels(expr$Language)=="en"] <- "English"
levels(expr$Language)[levels(expr$Language)=="de"] <- "German"
levels(expr$Language)[levels(expr$Language)=="fi"] <- "Finnish"
scatter_length <- ggplot(expr, aes(y=mean, x=length_bin, shape=Language, lty=Skill))

# plot of Tofd predicted by Length broken down by Language
# Figure 3 in manuscript
scatter_length + labs(x = "Length Bin", y = "Total Fixation Time")+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=position_dodge(width=0.1))+
  geom_point(position=position_dodge(width=0.1), size = 3, show.legend =T) + facet_grid( ~ Language)

### Plot of ToFD and word frequency for each language -- split by skill level
expr3 = summarySE(data = joint.fix.freq, measurevar = "ToFD", groupvars = c("Language", "Skill", "freq_bin"))
levels(expr3$Language)[levels(expr3$Language)=="en"] <- "English"
levels(expr3$Language)[levels(expr3$Language)=="de"] <- "German"
levels(expr3$Language)[levels(expr3$Language)=="fi"] <- "Finnish"
scatter_freq <- ggplot(expr3, aes(y=mean, x=freq_bin, shape=Language, lty=Skill))

# plot of Tofd x Frequency by Language
# Figure 4 in manuscript
scatter_freq + labs(x = "Frequency Bin", y = "Total Fixation Time")+
  geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), width=.2, position=position_dodge(width=0.1))+
  geom_point(position=position_dodge(width=0.1), size = 3, show.legend = T) + facet_grid( ~ Language)

####### Statistical Prediction for ToFD #######
de$lfreq = round(log(as.numeric(as.character(de$FREQcount))),2)
fi$lfreq = round(log(as.numeric(as.character(fi$FREQcount))),2)
en$lfreq = round(log(as.numeric(as.character(en$FREQcount))),2)

items = unique(en[, c("lfreq", "length")])
colnames(extremes) = c("ART", "SpellTest", "Lextale")

long_items = data.frame()
for (i in 1:8) {
  long_items = rbind(long_items, items)  
}
long_items$stam = paste(long_items$lfreq, long_items$length, sep = "_")
long_items[order(long_items$stam),] -> long_items
long_items$ART = extremes$ART
long_items$SpellTest = extremes$SpellTest
long_items$Lextale = extremes$Lextale
long_items -> dat
### Baseline models for German and Finnish 
# German Baseline Model
# Table A5
de.mod.cont = lmer(log(ToFD) ~ (length + lfreq) * (SpellTest + ART + Lextale) + 
                     (1 + length + lfreq| Subject) + (1 | TargetWord) + (1 | text), 
                   data = de)
summary(de.mod.cont)
cor(fitted(de.mod.cont), log(de$ToFD))^2

de.mod.cont1 = lmer(log(ToFD) ~ (length + lfreq) * (SpellTest + ART + Lextale) + 
                      (1 + length + lfreq| Subject) + (1 | TargetWord) + (1 | text), 
                    data = de, subset = abs(scale(resid(de.mod.cont)))<2.5)
summary(de.mod.cont1)

# Finnish Baseline Model
# Table A6
fi.mod.cont = lmer(log(ToFD) ~ (length + lfreq) * (SpellTest + ART + Lextale) +
                     (1 +length + lfreq| Subject) + (1 | TargetWord) + (1 | text), 
                   data = fi)
summary(fi.mod.cont)
cor(fitted(fi.mod.cont), log(fi$ToFD))^2

fi.mod.cont1 = lmer(log(ToFD) ~ (length + lfreq) * (SpellTest + ART + Lextale) + (1 + length + lfreq| Subject) + (1 | TargetWord) + (1 | text), 
                    data = fi, subset = abs(scale(resid(fi.mod.cont)))<2.5)
summary(fi.mod.cont1)
  
### Statistical Prediction Method
# Create the prediction model from the German baseline
predict(de.mod.cont1, newdata = dat, re.form = NA, allow.new.levels = T) -> pred1
dat -> dat.de
dat.de$fitted = pred1
# Create the prediction model from the Finnish baseline
predict(fi.mod.cont1, newdata = dat, re.form = NA, allow.new.levels = T) -> pred1
dat -> dat.fi
dat.fi$fitted = pred1

en$lToFD = log(en$ToFD)
de$lToFD = log(de$ToFD)
fi$lToFD = log(fi$ToFD)
dat.de$lToFD = dat.de$fitted
dat.fi$lToFD = dat.fi$fitted

#now we aggregate observed and predicted values by length and freq
dat.de$length_bin = cut(dat.de$length, breaks = c(1, 4, 5, 6,  8, 12), labels = c("3-4", "5", "6", "7-8", "9+"), right = T, include.lowest = T)
dat.fi$length_bin = cut(dat.fi$length, breaks = c(1, 4, 5, 6, 8, 12), labels = c("3-4", "5", "6", "7-8", "9+"), right = T, include.lowest = T)
de$length_bin = cut(de$length, breaks = c(1, 4, 5, 6, 8, 12), labels = c("3-4", "5", "6", "7-8", "9+"), right = T, include.lowest = T)
fi$length_bin = cut(fi$length, breaks = c(1, 4, 5, 6, 8, 12), labels = c("3-4", "5", "6", "7-8", "9+"), right = T, include.lowest = T)
en$length_bin = cut(en$length, breaks = c(1, 4, 5, 6, 8, 12), labels = c("3-4", "5", "6", "7-8", "9+"), right = T, include.lowest = T)

dat.de$freq_bin = cut(dat.de$lfreq, breaks = quantile(dat.de$lfreq), labels = c("low", "mid-low", "mid-high", "high"), right = T, include.lowest = T)
dat.fi$freq_bin = cut(dat.fi$lfreq, breaks = quantile(dat.fi$lfreq), labels = c("low", "mid-low", "mid-high", "high"), right = T, include.lowest = T)
de$freq_bin = cut(de$lfreq, breaks = quantile(dat.de$lfreq), labels = c("low", "mid-low", "mid-high", "high"), right = T, include.lowest = T)
fi$freq_bin = cut(fi$lfreq, breaks = quantile(dat.de$lfreq), labels = c("low", "mid-low", "mid-high", "high"), right = T, include.lowest = T)
en$freq_bin = cut(en$lfreq, breaks = quantile(dat.de$lfreq), labels = c("low", "mid-low", "mid-high", "high"), right = T, include.lowest = T)

#quantiles are aligned based on unique combinations of freq and length
table(dat.de$length_bin, dat.de$freq_bin)

dat.de = dat.de[!(dat.de$freq_bin == "high" & dat.de$length_bin == "9+"),]
dat.de = dat.de[!(dat.de$freq_bin == "low" & dat.de$length_bin == "3-4"),]
dat.fi = dat.fi[!(dat.fi$freq_bin == "high" & dat.fi$length_bin == "9+"),]
dat.fi = dat.fi[!(dat.fi$freq_bin == "low" & dat.fi$length_bin == "3-4"),]

de = de[!(de$freq_bin == "high" & de$length_bin == "9+"),]
de = de[!(de$freq_bin == "low" & de$length_bin == "3-4"),]
fi = fi[!(fi$freq_bin == "high" & fi$length_bin == "9+"),]
fi = fi[!(fi$freq_bin == "low" & fi$length_bin == "3-4"),]
en = en[!(en$freq_bin == "high" & en$length_bin == "9+"),]
en = en[!(en$freq_bin == "low" & en$length_bin == "3-4"),]

de.observed1 <- aggregate(lToFD ~ length_bin + freq_bin + Subject, data = de, FUN = median)
de.observed1a <- aggregate(lToFD ~ length_bin + freq_bin, data = de.observed1, FUN = median)
de.observed = rbind(de.observed1a)
de.observed$Type = "DE_OBS"

fi.observed1 <- aggregate(lToFD ~ length_bin + freq_bin + Subject, data = fi, FUN = median)
fi.observed1a <- aggregate(lToFD ~ length_bin + freq_bin, data = fi.observed1, FUN = median)
fi.observed = rbind(fi.observed1a)
fi.observed$Type = "FI_OBS"

en.observed1 <- aggregate(lToFD ~ length_bin + freq_bin + Subject, data = en, FUN = median)
en.observed1a <- aggregate(lToFD ~ length_bin + freq_bin, data = en.observed1, FUN = median)
en.observed1b <- aggregate(lToFD ~ length_bin + freq_bin, data = en.observed1, FUN = min)
en.observed1c <- aggregate(lToFD ~ length_bin + freq_bin, data = en.observed1, FUN = max)
en.observed = rbind(en.observed1a, en.observed1b, en.observed1c)
en.observed$Type = "EN_OBS"

de.predicted <- aggregate(lToFD ~ length_bin + freq_bin + ART + SpellTest + Lextale, data = dat.de, FUN = median)
de.predicted$Type = "DE"

fi.predicted <- aggregate(lToFD ~ length_bin + freq_bin + ART + SpellTest + Lextale, data = dat.fi, FUN = median)
fi.predicted$Type = "FI"

all = rbind(de.predicted[, c("length_bin",  "freq_bin", "lToFD", "Type")], 
            fi.predicted[, c("length_bin", "freq_bin", "lToFD", "Type")], 
            de.observed[, c("length_bin", "freq_bin", "lToFD", "Type")], 
            fi.observed[, c("length_bin", "freq_bin", "lToFD", "Type")], 
            en.observed[, c("length_bin", "freq_bin", "lToFD", "Type")])
table(all$Type)
all$Type = factor(all$Type, levels = c("DE_OBS", "DE", "EN_OBS", "FI", "FI_OBS"))
  
### Plot showing all predictions of length x frequency combinations of hypothetical Finnish and German speakers
all$factor = factor(paste(all$length_bin, all$freq_bin, sep = "_"))
# Formatting to get plotting points
all_obs <- all %>% filter(Type %in% c("DE_OBS", "FI_OBS", "EN_OBS")) %>% droplevels()
all_obs <- all_obs %>% 
  group_by(factor) %>% 
  mutate(todf_type = ifelse(lToFD == min(lToFD), "min",
                            ifelse(lToFD == max(lToFD), "max", lToFD))) %>% 
  filter(!(todf_type %in% c("min", "max"))) %>% droplevels()  
all_pred_fi <- all %>% filter(Type == "FI") %>% droplevels()
all_pred_de <- all %>% filter(Type == "DE") %>% droplevels()
all_obs$factor <- forcats::fct_relevel(all_obs$factor,
                                         "3-4_mid-low",
                                         "3-4_mid-high",
                                         "3-4_high",
                                         "5_low", 
                                         "5_mid-low",
                                         "5_mid-high",
                                         "5_high",
                                         "6_low",
                                         "6_mid-low",
                                         "6_mid-high",
                                         "6_high",
                                         "7-8_low",
                                         "7-8_mid-low",
                                         "7-8_mid-high",
                                         "7-8_high",
                                         "9+_low",
                                         "9+_mid-low",
                                         "9+_mid-high")
cols <- c("EN_OBS" = "black", "FI_OBS" = "black", "FI" = "grey40")

# Figure 5: Finnish observed, English observed, Finnish predicted
all_obs %>% 
  filter(Type %in% c("EN_OBS", "FI_OBS")) %>% 
  ggplot(aes(factor, (exp(lToFD)), color = Type, group = Type, shape = Type)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(17, 1, 1), guide = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #geom_line() +
  geom_point(data = all_pred_fi, alpha = 0.5, shape = 20, size = 3) +
  labs(x = "Length and frequency bin", y = "Total fixation duration") + 
  scale_colour_manual(values = cols) + 
  #scale_color_manual(values = c("grey40", "grey80", "black")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12,
                                   angle = 90),
        axis.text.y = element_text(size = 12),
        legend.position ="bottom",
        legend.box = "horizontal")

# Figure 6: German observed, English observed, German predicted
cols2 <- c("EN_OBS" = "black", "DE_OBS" = "black", "DE" = "grey40")
all_obs %>% 
  filter(Type %in% c("EN_OBS", "DE_OBS")) %>% 
  ggplot(aes(factor, (exp(lToFD)), color = Type, group = Type, shape = Type)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(17, 5, 5), guide = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #geom_line() +
  geom_point(data = all_pred_de, alpha = 0.5, shape = 20, size = 3) +
  labs(x = "Length and frequency bin", y = "Total fixation duration") + 
  scale_colour_manual(values = cols2) +
  #scale_color_manual(values = c("grey40", "grey80", "black")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 12,
                                   angle = 90),
        axis.text.y = element_text(size = 12),
        legend.position ="bottom",
        legend.box = "horizontal")
  
####### Cognates #######
# read in file containing the list of cognates in each language
cognates = read.csv("words_cognate_DeFi.csv")
# Data formatting and cleaning
cognates = dplyr::select(cognates, "TargetWord","Cognate_De", "Cognate_Fi")
cognates$TargetWord = gsub("\\.", "", cognates$TargetWord)
# remove stopwords
cognates = cognates[!cognates$TargetWord %in% stopwords,] 
# creating a list of just the unique words when adding the cognates to joint.fix.cog
cognates = unique(cognates[,c("TargetWord","Cognate_De", "Cognate_Fi")]) 
cognates = cognates[-93,] # removing an incorrectly identified cognate
# joining the cognate data with joint.fix
joint.fix.cog = merge(joint.fix.freq, cognates, by.x = "TargetWord", by.y = "TargetWord")
# 5 words found in joint.fix but not in the list of cognates are removed
setdiff(joint.fix$TargetWord, cognates$TargetWord)
setdiff(cognates$TargetWord, joint.fix$TargetWord)
length(cognates$TargetWord)

# Number of cognates in each language
table(cognates$Cognate_De)
# 156 / 316
table(cognates$Cognate_Fi)
# 44 / 316

### No effect for German
lm_cog_De <- lmer(log(ToFD) ~ Cognate_De + Skill  + log(SUBTLEX) + length_bin + (1 | TargetWord) + 
                    (1 | Subject), data = joint.fix.cog[joint.fix.cog$Language == "de",])
summary(lm_cog_De)

### No effect for Finnish
lm_cog_fi <- lmer(log(ToFD) ~ Cognate_Fi * Skill + log(SUBTLEX)  + length_bin +  (1 | TargetWord) + (1 | Subject), data = joint.fix.cog[joint.fix.cog$Language == "fi",])
summary(lm_cog_fi)

####### Component Skills #######
# pick out all of the predicted hypothetical speakers that did the best
# For German, these are all those that exceeded the highest observed English scores
# For Finnish, these are all those that exceeded the 40% percentil of observed English scores
u = levels(all$factor)
count.de = count.fi  = 0
lis.de = lis.fi = list()
for (i in 1:length(u)) {
  tmp = all[all$factor == u[i],]
  w = which(tmp[tmp$Type == "DE",]$lToFD < quantile(tmp[tmp$Type == "EN_OBS",]$lToFD, 0.1))
  w.fi = which(tmp[tmp$Type == "FI",]$lToFD <= quantile(tmp[tmp$Type == "EN_OBS",]$lToFD, 0.6)) #median is too hard. Trying 40%
  if (length(w)>0) {
    count.de = count.de + 1
    tmp1 = tmp[tmp$Type == "DE",]
    print(tmp1[w,])
    lis.de[[count.de]] = tmp1[w,]
  }
  
  if (length(w.fi)>0) {
    count.fi = count.fi + 1
    tmp1 = tmp[tmp$Type == "FI",]
    print(tmp1[w.fi,])
    lis.fi[[count.fi]] = tmp1[w.fi,]
  }
}

fast.de = data.frame(do.call(rbind, lis.de))
fast.fi = data.frame(do.call(rbind, lis.fi))

fast.de <- merge(fast.de, dat.de)
fast.fi <- merge(fast.fi, dat.fi)

### Top predicted combinations for German
# Table 7
fast.de
### Top predicted combinations for Finnish
# Table 8
fast.fi