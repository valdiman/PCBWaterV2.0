## Pre water PCB concentrations data analysis

# Install packages
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("robustbase")
install.packages("dplyr")
install.packages("tibble")
install.packages("Matrix")
install.packages("lme4")
install.packages("MuMIn")
install.packages("lmerTest")

# Load libraries
library(ggplot2)
library(scales) # function trans_breaks
#library(gridExtra)
#library(tidyverse)
library(stringr) # str_detect
library(robustbase) #function colMedians
library(dplyr) # %>%
library(tibble) # function add a column
library(lme4) # perform lme
library(MuMIn) # get Rs from lme
library(lmerTest) # get the p-value from lme

# Read data ---------------------------------------------------------------

# Data in pg/L
wdc.0 <- read.csv("WaterDataCongenerAroclor062322.csv")

# Data analysis -----------------------------------------------------------

# Look at potential distributions
# Needs to add a individual number for each site name
# Create a site number
site.numb <- wdc.0$SiteName %>% as.factor() %>% as.numeric
wdc.site <- add_column(wdc.0,
                       site.numb, .after = "AroclorCongener")
wdc.site$SampleDate <- as.Date(wdc.site$SampleDate,
                               format = "%m/%d/%y")
# Calculate total PCB per sample
tpcb <- rowSums(wdc.site[, c(13:116)], na.rm = T)
time.day <- data.frame(as.Date(wdc.site$SampleDate) - min(as.Date(wdc.site$SampleDate)))
# Generate data.frame for analysis and plots
tpcb <- cbind(data.frame(time.day), as.matrix(tpcb),
                  wdc.site$SampleDate)
colnames(tpcb) <- c("time", "tPCB", "date")

# Calculate total log PCB per sample
# Remove metadata
wdc.site.1 <- subset(wdc.site, select = -c(ID:site.numb))
# Remove Aroclor data
wdc.site.1 <- subset(wdc.site.1, select = -c(A1016:A1260))
log.pcb <- log10(wdc.site.1)
t.log.pcb <- rowSums(log.pcb, na.rm = T)
# Generate data.frame for analysis and plots
log.tpcb <- cbind(data.frame(time.day), as.matrix(t.log.pcb),
              wdc.site$SampleDate)
colnames(log.tpcb) <- c("time", "tPCB", "date")

# Histograms
hist(tpcb$tPCB)
hist(log10(tpcb$tPCB))
hist(log.tpcb$tPCB)

# Regressions
# (1) Total PCB, tpcb
# (1.1) Perform linear regression (lr)
# + 1
lr.tpcb <- lm(log10(tPCB + 1) ~ time, data = tpcb)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (1.2) Perform Linear Mixed-Effects Model (LME)
# Site number code
site <- wdc.site$site.numb
time <- tpcb$time
# + 1
lmem.tpcb <- lmer(log10(tpcb$tPCB + 1) ~ 1 + time + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.tpcb)
# Look at residuals
res <- resid(lmem.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res, main = "log10(C + 1)")
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2) Sum of log10 individual PCBs, log.tpcb
# (2.1) Perform linear regression (lr)
# + 1
lr.tpcb <- lm(log10(tPCB + 1) ~ time, data = log.tpcb)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# (2.2) Perform Linear Mixed-Effects Model (LME)
# Site number code
site <- wdc.site$site.numb
time <- tpcb$time
# + 1
lmem.tpcb <- lmer(log10(log.tpcb$tPCB + 1) ~ 1 + time + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.tpcb)
# Look at residuals
res <- resid(lmem.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res, main = "log10(C + 1)")
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')

# Analysis per site -------------------------------------------------------
# Per site

Blu.R <- wdc.0[str_detect(wdc.0$SiteName, 'BlueRiver'),]
Fox.R <- wdc.0[str_detect(wdc.0$SiteName, 'FoxRiver'),]
Hud.R <- wdc.0[str_detect(wdc.0$SiteName, 'HudsonRiver'),]
Hou.R <- wdc.0[str_detect(wdc.0$SiteName, 'HousatonicRiver'),]
Kal.R <- wdc.0[str_detect(wdc.0$SiteName, 'KalamazooRiver'),]





# Separated by dates
tpcb.1995 <- tpcb[tpcb$date <= 1995, ]
tpcb.2000 <- tpcb[tpcb$date > 1995 & tpcb$date <= 2000, ]
tpcb.2005 <- tpcb[tpcb$date > 2000 & tpcb$date <= 2005, ]
tpcb.2010 <- tpcb[tpcb$date > 2005 & tpcb$date <= 2010, ]
tpcb.2015 <- tpcb[tpcb$date > 2010 & tpcb$date <= 2015, ]
tpcb.2020 <- tpcb[tpcb$date > 2015 & tpcb$date <= 2020, ]
tpcb.2000.1 <- tpcb[tpcb$date <= 2000, ]
tpcb.2020.1 <- tpcb[tpcb$date > 2000 & tpcb$date <= 2020, ]
tpcb.2005.1 <- tpcb[tpcb$date <= 2005, ]
tpcb.2020.2 <- tpcb[tpcb$date > 2005 & tpcb$date <= 2020, ]



hist(tpcb.1995$tPCB)
hist(log10(tpcb.1995$tPCB))
# Create Q-Q plot for residuals
qqnorm(log10(tpcb.1995$tPCB + 1))
# Add a straight diagonal line to the plot
qqline(log10(tpcb$tPCB + 1))

hist(tpcb.2000$tPCB)
hist(log10(tpcb.2000$tPCB))
# Create Q-Q plot for residuals
qqnorm(log10(tpcb.2000$tPCB + 1))
# Add a straight diagonal line to the plot
qqline(log10(tpcb.2000$tPCB + 1))

hist(tpcb.2005$tPCB)
hist(log10(tpcb.2005$tPCB))
# Create Q-Q plot for residuals
qqnorm(log10(tpcb.2005$tPCB + 1))
# Add a straight diagonal line to the plot
qqline(log10(tpcb.2005$tPCB + 1))

hist(tpcb.2010$tPCB)
hist(log10(tpcb.2010$tPCB))
# Create Q-Q plot for residuals
qqnorm(log10(tpcb.2010$tPCB + 1))
# Add a straight diagonal line to the plot
qqline(log10(tpcb.2010$tPCB + 1))

hist(tpcb.2015$tPCB)
hist(log10(tpcb.2015$tPCB))
# Create Q-Q plot for residuals
qqnorm(log10(tpcb.2015$tPCB + 1))
# Add a straight diagonal line to the plot
qqline(log10(tpcb.2015$tPCB + 1))

hist(tpcb.2020$tPCB)
hist(log10(tpcb.2020$tPCB))
# Create Q-Q plot for residuals
qqnorm(log10(tpcb.2020$tPCB + 1))
# Add a straight diagonal line to the plot
qqline(log10(tpcb.2020$tPCB + 1))

# Perform linear regression (lr)
# + 1
lr.tpcb <- lm(log10(tPCB + 1) ~ time, data = tpcb.2005.1)
# See results
summary(lr.tpcb)
# Look at residuals
res <- resid(lr.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res)
# Add a straight diagonal line to the plot
qqline(res)
# Shapiro test
shapiro.test(res)
# One-sample Kolmogorov-Smirnov test
ks.test(res, 'pnorm')
