# Required libraries
library(readxl)
library(tidyverse)
library(hrbrthemes)
library(lme4)
library(ggeffects)
library(sjPlot)

setwd('~/Dropbox/Projects/GBM/Codes/DIPG/DIPG_RT_drug')

# Import data
dat <- read_excel("DIPG6_BLI_dat.xlsx")

###############################################################################

dat$Tumor_vol_log = log(dat$Tumour_Vol)

head(dat)

ggplot(dat, aes(time, Tumor_vol_log, MouseID = Mouse_ID, col = Cohort)) +
  geom_path(size = 1) +
  scale_color_viridis_d() +
  # scale_colour_manual(values=c("#404788FF","#2D708EFF","#33638DFF","20A387FF")) +
  # geom_point() +
  theme_classic() +
  labs(x = "Time (days)", y = "Radiance (log) (p/s)")

