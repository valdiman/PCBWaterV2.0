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
wdc.0 <- read.csv("WaterDataCongenerAroclor072122.csv")

# Prepare data -----------------------------------------------------------

# Look at potential distributions
# (i) Calculate total PCB per sample
tpcb <- rowSums(wdc.0[, c(12:115)], na.rm = T)
# Change date format
wdc.0$SampleDate <- as.Date(wdc.0$SampleDate, format = "%m/%d/%y")
# Calculate sampling time
time.day <- data.frame(as.Date(wdc.0$SampleDate) - min(as.Date(wdc.0$SampleDate)))
# Create data frame
tpcb <- cbind(wdc.0$SiteName, wdc.0$SampleDate,
              as.matrix(tpcb), data.frame(time.day), wdc.0$AroclorCongener)
# Add names to columns
colnames(tpcb) <- c("Site", "date", "tPCB", "time", "AroclorCongener")

# (ii) Calculate total log PCB per sample
# Remove metadata
log.pcb <- subset(wdc.0, select = -c(ID:AroclorCongener))
# Remove Aroclor data
log.pcb <- subset(log.pcb, select = -c(A1016:A1260))
# Log 10 individual PCBs 
log.pcb <- log10(log.pcb)
# Sum individual log 10 PCBs
t.log.pcb <- rowSums(log.pcb, na.rm = T)
# Replace -inf to NA
t.log.pcb[is.infinite(t.log.pcb)] = NA
# Generate data.frame for analysis and plots
log.tpcb <- cbind(wdc.0$SiteName, wdc.0$SampleDate,
                  as.matrix(t.log.pcb), data.frame(time.day),
                  wdc.0$AroclorCongener)
colnames(log.tpcb) <- c("Site", "date", "logtPCB", "time", "AroclorCongener")

# Histograms --------------------------------------------------------------

# (i)
hist(tpcb$tPCB)
hist(log10(tpcb$tPCB))
hist(tpcb$tPCB[tpcb$tPCB < 1000])
hist(log10(tpcb$tPCB[tpcb$tPCB < 10^3]))
# (ii)
hist(log.tpcb$logtPCB)
hist(log10(log.tpcb$logtPCB))
hist(log10(1 + 122.22764 + log.tpcb$logtPCB))

# Regressions -------------------------------------------------------------

# (i) Total PCB, tpcb
# (i.1) Perform linear regression (lr)
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

# (i.2) Perform Linear Mixed-Effects Model (LME)
# Site number code
# Needs to add a individual number for each site name
# Create a site number
site.numb <- wdc.0$SiteName %>% as.factor() %>% as.numeric
wdc.site <- add_column(wdc.0, site.numb, .after = "AroclorCongener")
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

# (ii) Sum of log10 individual PCBs, log.tpcb
# (ii.1) Perform linear regression (lr)
lr.tpcb <- lm(logtPCB ~ time, data = log.tpcb)
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

# (ii.2) Perform Linear Mixed-Effects Model (LME)
lmem.tpcb <- lmer(log10(1 + 122.22764 + log.tpcb$logtPCB) ~ 1 + time + (1|site),
                  REML = FALSE,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nRE="ignore"))

# See results
summary(lmem.tpcb)
# Look at residuals
res <- resid(lmem.tpcb) # get list of residuals
# Create Q-Q plot for residuals
qqnorm(res, main = "C")
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

# Select site to work
wdc.s <- Kal.R
# (i) Calculate total PCB per sample
tpcb <- rowSums(wdc.s[, c(12:115)], na.rm = T)
# Change date format
wdc.s$SampleDate <- as.Date(wdc.s$SampleDate, format = "%m/%d/%y")
# Calculate difference between initial time and rest
time.day <- data.frame(as.Date(wdc.s$SampleDate) - min(as.Date(wdc.s$SampleDate)))
# Generate data.frame for analysis and plots
tpcb.s <- cbind(wdc.s$SampleDate, as.matrix(tpcb), data.frame(time.day))
colnames(tpcb.s) <- c("date", "tPCB", "time")

# (ii) Calculate total log PCB per sample
# Remove metadata
wdc.s.1 <- subset(wdc.s, select = -c(ID:AroclorCongener))
# Remove Aroclor data
wdc.s.1 <- subset(wdc.s.1, select = -c(A1016:A1260))
# Log 10 individual PCBs 
log.pcb.s <- log10(wdc.s.1)
# Sum individual log 10 PCBs
t.log.pcb.s <- rowSums(log.pcb.s, na.rm = T)
# Replace -inf to NA
t.log.pcb.s[is.infinite(t.log.pcb.s)] = NA
# Generate data.frame for analysis and plots
log.tpcb.s <- cbind(wdc.s$SampleDate, as.matrix(t.log.pcb.s),
                    data.frame(time.day))
colnames(log.tpcb.s) <- c("date", "logtPCB", "time")

# Histograms
hist(tpcb.s$tPCB)
hist(log10(tpcb.s$tPCB))
hist(log.tpcb.s$logtPCB)
hist(log10(log.tpcb.s$logtPCB))

# Regressions
# (1) Total PCB, tpcb
# (1.1) Perform linear regression (lr)
# + 1
lr.tpcb <- lm(log10(tPCB + 1) ~ time, data = tpcb.s)
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

# (2) Sum of log10 individual PCBs, log.tpcb
# (2.1) Perform linear regression (lr)
lr.tpcb <- lm(logtPCB ~ time, data = log.tpcb.s)
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

# Data per date -----------------------------------------------------------

# Separated data per dates
# (i) PCBs
# Create column with just year
tpcb$year <- format(as.Date(tpcb$date, format="%d/%m/%Y"),"%Y")
# Years
tpcb.1995 <- tpcb[tpcb$year <= 1995, ]
tpcb.2000 <- tpcb[tpcb$year > 1995 & tpcb$year <= 2000, ]
tpcb.2005 <- tpcb[tpcb$year > 2000 & tpcb$year <= 2005, ]
tpcb.2010 <- tpcb[tpcb$year > 2005 & tpcb$year <= 2010, ]
tpcb.2015 <- tpcb[tpcb$year > 2010 & tpcb$year <= 2015, ]
tpcb.2020 <- tpcb[tpcb$year > 2015 & tpcb$year <= 2020, ]
tpcb.2000.1 <- tpcb[tpcb$year <= 2000, ]
tpcb.2020.1 <- tpcb[tpcb$year > 2000 & tpcb$year <= 2020, ]
tpcb.2005.1 <- tpcb[tpcb$year <= 2005, ]
tpcb.2020.2 <- tpcb[tpcb$year > 2005 & tpcb$year <= 2020, ]

# Histograms
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

# (ii) log10 PCBs
# Create column with just year
tpcb$year <- format(as.Date(log.tpcb$date, format="%d/%m/%Y"),"%Y")
# Years
tpcb.1995 <- tpcb[tpcb$year <= 1995, ]
tpcb.2000 <- tpcb[tpcb$year > 1995 & tpcb$year <= 2000, ]
tpcb.2005 <- tpcb[tpcb$year > 2000 & tpcb$year <= 2005, ]
tpcb.2010 <- tpcb[tpcb$year > 2005 & tpcb$year <= 2010, ]
tpcb.2015 <- tpcb[tpcb$year > 2010 & tpcb$year <= 2015, ]
tpcb.2020 <- tpcb[tpcb$year > 2015 & tpcb$year <= 2020, ]
tpcb.2000.1 <- tpcb[tpcb$year <= 2000, ]
tpcb.2020.1 <- tpcb[tpcb$year > 2000 & tpcb$year <= 2020, ]
tpcb.2005.1 <- tpcb[tpcb$year <= 2005, ]
tpcb.2020.2 <- tpcb[tpcb$year > 2005 & tpcb$year <= 2020, ]

# Histograms
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

