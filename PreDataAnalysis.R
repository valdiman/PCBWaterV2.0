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

# Look at potential distributions
# Needs to add a individual number for each site name
# Create a site number
site.numb <- wdc.0$SiteName %>% as.factor() %>% as.numeric
# Include site number after AroclorCongener column
wdc.site <- add_column(wdc.0,
                       site.numb, .after = "AroclorCongener")
# Format date
wdc.site$SampleDate <- as.Date(wdc.site$SampleDate,
                               format = "%m/%d/%y")
# Calculate total PCB per sample
tpcb <- rowSums(wdc.site[, c(13:116)], na.rm = T)
# Create sampling days from the first sample date
time.day <- data.frame(as.Date(wdc.site$SampleDate) - min(as.Date(wdc.site$SampleDate)))
# Generate data.frame w/ time.day, tpcb and modify date
tpcb <- cbind(data.frame(time.day), as.matrix(tpcb),
                  wdc.site$SampleDate)
colnames(tpcb) <- c("time", "tPCB", "date")



