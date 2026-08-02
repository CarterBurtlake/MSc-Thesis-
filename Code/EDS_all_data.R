# RLS + DFO + Storm coast + Christine Hansen + Jane Watson data
# Carter Burtlake 
# August 2026

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

##I am going to work in chronological order too pull all the data I will use in this project into this script, and then join them all into a unified df

#Jane Watson data ------------------------------------------------------------------------
JW_data <- read.csv("data-raw/AbaloneJaneWatsonCombined89_22.csv")

#lets review it
View(JW_data)
#all sites accounted for (84)
#TO DO: 1. column titles 2.site names (remove numbers) 3. sizes 4. sample size (remove n=)
# 5. update density relative to 1m^2 
JW_data_clean <- clean_names(JW_data) %>%
  mutate( #tidy site names
    code = as.integer(str_extract(site, "^\\d+")),
    site = case_when(
      code == 1 ~ "ship_islands",
      code == 2 ~ "ed_king",
      code == 3 ~ "w_ed_king",
      code == 4 ~ "seppings",
      code == 5 ~ "cape_beale",
      code == 6 ~ "lawton_point",
      code == 7 ~ "whittlestone",
      code == 8 ~ "execution_rock",
      code == 9 ~ "self_point",
      code == 10 ~ "helby_island",
      code == 11 ~ "scotts_bay",
      code == 12 ~ "aguilar_point",
      code == 13 ~ "cia_rock",
      code == 14 ~ "taylor_rock",
      code == 15 ~ "kirby_point",
      code == 16 ~ "village_bay",
      code == 17 ~ "blowhole",
      code == 18 ~ "prasiola",
      code == 19 ~ "wizard",
      code == 20 ~ "grappler",
      code == 21 ~ "ed_king_?",
      TRUE ~ site
    )
  ) %>%
  select(-code)%>%
  mutate(density_sample_size = str_remove(density_sample_size, "^n="))%>% #remove n= from rows in density sample size
  mutate(mean_density_1m = mean/0.49)%>% #okay so quadrats are 0.7m x 0.7m = 0.49m^2, so, to get abalone mean in quadrats across 1m2 we must do mean abalone in 0.49m2 / 0.49m2 = abalone per m2
  mutate(depth_m = depth_ft * 0.3048)%>%
  mutate(date = mdy(date), year = year(date)) #changes all the date formats using lubridate so that the second call of "year" can be applied to make a new column with just year

  
#we are going to leave the size information for JW right now as really dont have much of that...
#no depth for 1994 - perhaps it was assumed to be the same as previous survey?
#also note depth makes major chage throughout program

View(JW_data_clean)
#readr::write_csv(JW_data_clean, file = "data-processed/JW_data_1988_2022.csv") #to save the clean data

#Christine Hansen data -------------------------------------------------------------------
#in the future I hope to have this raw data...
#for now, to get this data we have to use meta digitize 
#CH_data_shallow <- metaDigitise("/Users/carterburtlake/Documents/RProjects/MScChapter2/figures")
#load in shallow data

#CH_data_deep <- metaDigitise("/Users/carterburtlake/Documents/RProjects/MScChapter2/figures")
#load in deep data

#readr::write_csv(CH_data_deep, file = "data-raw/CH_data_2008.csv") #to save the raw data

#TO DO: rename columns to match JW dataset then save as data
#need site, date, depth, add zero counts, mean_density_1m

CH_data_clean <- CH_data_deep%>%
  mutate(study = "ChristineHansen")%>% #too keep df's consistent
  mutate(year = "2008")%>% #add year data from paper
  mutate(depth_m = case_when(str_detect(group_id, "deep") ~"7.5", TRUE ~ "2.5")) %>%
  #from the paper we know that the shallow sites are from 0-5m and deep sites are 5-10 so just take the middle of both
  #revist this decision for more complicated analyses 
  mutate(reproductive_status = case_when(str_detect(group_id, "_i") ~ "immature", TRUE ~ "reproductive")) %>%
  mutate(site = case_when(str_detect(group_id, "ellis") ~ "ellis",
                          str_detect(group_id, "helby") ~ "helby",
                          TRUE ~ "scotts_bay")) %>% #give them independent site names outside of depth and reproductive status
  add_row(study = "ChristineHansen", mean = 0, se = 0, year = "2008", depth_m = "7.5", reproductive_status = "reproductive", site = "helby")%>%
  add_row(study = "ChristineHansen", mean = 0, se = 0, year = "2008", depth_m = "2.5", reproductive_status = "reproductive", site = "helby")%>%
  add_row(study = "ChristineHansen", mean = 0, se = 0, year = "2008", depth_m = "7.5", reproductive_status = "reproductive", site = "goby_town")%>%
  add_row(study = "ChristineHansen", mean = 0, se = 0, year = "2008", depth_m = "7.5", reproductive_status = "immature", site = "goby_town")%>%
  add_row(study = "ChristineHansen", mean = 0, se = 0, year = "2008", depth_m = "2.5", reproductive_status = "reproductive", site = "goby_town")%>%
  add_row(study = "ChristineHansen", mean = 0, se = 0, year = "2008", depth_m = "2.5", reproductive_status = "immature", site = "goby_town") #had to add all the goby town data where they searched for abalone but found none during their surveys. Same is also true from "reproductive" abalone
#note that year, and depth are characters... we will probably need them as numeric and dates at some point so lets change that now
is.character(CH_data_clean$year) #TRUE
is.character(CH_data_clean$depth_m) #TRUE


CH_data_mean <- CH_data_clean%>%
  group_by(site, year)%>%
  summarise(site_mean = sum(mean, na.rm = TRUE))
#add the means across the whole site regardless of depth or size class / reproductive status 
#this is processed data and should be saved as such

#readr::write_csv(CH_data_mean, file = "data-processed/CH_data_2008.csv") #to save the processed data

#DFO kelp data ------------------------------------------------------------------------
DFO_kelp_data <- read.csv("data-raw/2021-2025_DFO_Kelp_Abalone.csv")
