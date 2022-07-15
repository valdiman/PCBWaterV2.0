## Pre water PCB concentrations data analysis

# Install packages
install.packages("ggpubr")
install.packages("ggpmisc")
install.packages("tidyverse")
install.packages("reshape2")
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
#library(ggpubr)
#library(ggpmisc)
library(scales) # function trans_breaks
#library(gridExtra)
#library(tidyverse)
#library(reshape2)
library(stringr)
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

# Needs to add a individual number for each site name
# to perform regression analysis, linear and Linear Mixed-Effects Model (LME)
site.numb <- wdc.0$SiteName %>% as.factor() %>% as.numeric
wdc.site <- add_column(wdc.0,
                       site.numb, .after = "AroclorCongener")
wdc.site$SampleDate <- as.Date(wdc.site$SampleDate,
                               format = "%m/%d/%y")
# Calculate total PCB per sample
tpcb.tmp <- rowSums(wdc.site[, c(13:116)], na.rm = T)
time.day <- data.frame(as.Date(wdc.site$SampleDate) - min(as.Date(wdc.site$SampleDate)))
# Generate data.frame for analysis and plots
tpcb.tmp <- cbind(data.frame(time.day), as.matrix(tpcb.tmp),
                  wdc.site$SampleDate)
colnames(tpcb.tmp) <- c("time", "tPCB", "date")