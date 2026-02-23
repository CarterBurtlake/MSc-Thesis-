# RLS abalone data
# Carter Burtlake 
# February 2026

# Load packages and data -----
# Load packages
library(readr)
library(tidyr)
library(dplyr) # for general data wrangling
library(lubridate) # for looking at dates

library(glmmTMB) # The swiss army knife of modeling packages
library(DHARMa) # inspect model residuals/check assumptions
library(ggeffects) # for extracting predictions and running post hoc tests

library(ggplot2) # Plotting data
library(patchwork) # Arrange multiple plots together
library(visreg) # plot model predictions
library(viridisLite) # colour palettes

library(ggspatial)
library(sf)

# Load and manipulate data -----------------------------------------------------
abalone <- read_csv("Data/abalone_RLS.csv") # this is the averaged abalone data from RLS

