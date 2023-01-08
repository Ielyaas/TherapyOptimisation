# Required libraries
library(readxl)
library(tidyverse)
library(hrbrthemes)
library(lme4)
library(ggeffects)
library(sjPlot)

setwd('~/Dropbox/Projects/GBM/Codes/DIPG/DIPG_RT_drug')

# Import data
dat <- read_excel("BLI_dat_processed_150822.xlsx",sheet = "DIPG7_combine")

###############################################################################

dat$Tumor_vol_log = log(dat$Tumour_Vol)

head(dat)

ggplot(dat, aes(time, Tumor_vol_log, MouseID = Mouse_ID, col = Cohort, linetype = as.factor(Experiment))) +
  geom_path() +
  # geom_point() +
  facet_wrap(~factor(Cohort)) +
  labs(x = "Time (days)", y = "Radiance (log) (p/s)")

