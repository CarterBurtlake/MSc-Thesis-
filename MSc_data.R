#MSc data
#Carter Burtlake
#October 2026

# Load packages and data -----
# Load packages
library(readr)
library(tidyr)
library(dplyr) # for general data wrangling
library(lubridate) # for looking at dates
library(janitor) #to clean data frames
library(stringr) #to clean data frames

library(glmmTMB) # The swiss army knife of modeling packages
library(DHARMa) # inspect model residuals/check assumptions
library(ggeffects) # for extracting predictions and running post hoc tests

library(ggplot2) # Plotting data
library(patchwork) # Arrange multiple plots together
library(visreg) # plot model predictions
library(viridisLite) # colour palettes
library(metaDigitise) #to pull data from figures without raw values

library(ggspatial)
library(sf)

# Load DFO kelp data ----
DFO_kelp_data <- read.csv("data-raw/2021-2025_DFO_Kelp_Abalone.csv")
DFO_site_data <- read.csv("data-raw/DFO_Kelp_Sites.csv")

#hmm this data set looks different... I need to get year and site into this data frame...
#From the methods I know the area covered on each survey, note that some comments or notes might describe missing transect lines or quadrats

#Filter for just abalone
DFO_kelp_data_clean <- DFO_kelp_data %>%
  filter(SpeciesCode == "A") %>%
  uncount(Count) %>% #change rows that only contain counts and not measures (e.g., measure = N/A, Count = 7 -> 7 new rows with same parent ID and no measurement info in each row)
  mutate(Count = 1)

#Add the site data using the unique ID for each row
DFO_abalone <- DFO_site_data %>%
  left_join(DFO_kelp_data_clean %>%
              select(HKey, Quadrat, GaugeDepthFeet, Number, SpeciesCode, Measurement, Count, Comments),
            by = "HKey")

#Next we need to add a total for each transect (WHICH I MUST CONFIRM ALL USED 3 QUADRATS TO GET THAT VALUE) for each year, add them up across the whole site, and then divide by the appropriate area of size surveyed to get density for a year. 
