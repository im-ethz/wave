# Packages
library(tidyverse)
library(readr)
library(readxl)

library(abind)
library(magrittr)
library(sjmisc)
library(RcppRoll)
library(zoo)
library(lubridate)
library(hms)

library(glmnet)
library(brms)
options(mc.cores = parallel::detectCores())
library(tidybayes)
library(loo)

library(survival)
library(survminer)

library(beepr)

# Plotting
library(ggplot2)
library(plotly)
library(scales)
library(ggdist)
library(ggrepel)
library(ggsci)
library(cowplot)
library(bayesplot)
library(gginnards)

discharge_color <- "#4DBBD5FF"
icu_color <- "#E64B35FF"

white_plot <- ggplot(data.frame()) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "white"))

# Fonts
library(extrafont)

# Note: an older version of Rttf2pt1 is currently required for the fonts to work properly
# library(remotes)
# remotes::install_version("Rttf2pt1", version = "1.3.8")

font_import(pattern = "times", prompt = F)
loadfonts(device = "win")
loadfonts(device = "postscript")
loadfonts(device = "pdf")
Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.55.0/bin/gswin64c.exe") # ghostscript
