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

# DFO kelp data ----
#load data
DFO_kelp_data <- read.csv("data-raw/2021-2025_DFO_Kelp_Abalone.csv")
DFO_site_data <- read.csv("data-raw/DFO_Kelp_Sites.csv")

#Filter for just abalone
DFO_kelp_data_clean <- DFO_kelp_data %>%
  filter(SpeciesCode == "A") %>% #this is the species code shared in the DFO database
  uncount(Count) %>% #change rows that only contain counts and not measures (e.g., measure = N/A, Count = 7 -> 7 new rows with same parent ID and no measurement info in each row)
  mutate(Count = 1)

#Add the site data using the unique ID for each row
DFO_abalone <- DFO_site_data %>%
  left_join(DFO_kelp_data_clean %>%
              select(HKey, Quadrat, GaugeDepthFeet, Number, SpeciesCode, Measurement, Count, Comments),
            by = "HKey")

#Next we need to add a total for each transect (WHICH I MUST CONFIRM ALL USED 3 QUADRATS TO GET THAT VALUE) for each year, add them up across the whole site, and then divide by the appropriate area of size surveyed to get density for a year (mean). 

DFO_abalone_density <- DFO_abalone %>%
  group_by(Year, SiteName) %>%
  summarise(total_abalone = sum(Count, na.rm = TRUE),
            n_transects = n_distinct(Transect[!is.na(Transect)]), #to capture the number of transects that go into a site for a particular year
            .groups = "drop") %>%
  mutate(Area = n_transects * 3) %>% #take the number of transects at that site and * them by the number of quadrats used along a transect = 3 *COFNRIM THIS WITH SANDIE*
  mutate(Density = total_abalone / Area) #calc density

#lets get an average for each site across all years 2021-2025
DFO_average_density <- DFO_abalone_density %>%
  filter(!SiteName %in% c("Taylor Island 2021", "Flemming 112 2021", "Flemming 114")) %>% #remove sites that were only done in 2021
  # CHECK: Confirm we should remove or keep these? I think we should remove in this specific case as the sites were seen as unsuitable
  group_by(SiteName) %>%
  summarise(mean_density = mean(Density), #calculate density average across all years of study
            number_years = n_distinct(Year), #tell us how many years go into the average calculation
            .groups = "drop")



# RLS data ----
#load data
RLS_data <- read.csv("data-raw/RLS_2025.csv") 

#filter for sites over years that don't have abalone 
RLS_data_year <- RLS_data %>%
  mutate(survey_date = ymd(survey_date)) %>%
  mutate(year = year(survey_date)) #first we make a year column

#attempt to clean this for surveys in a specific year, at a site, at a depth that do or dont have abalone
# All unique site-year-depth combinations
all_combos <- RLS_data_year %>%
  distinct(site_code, year, depth)

# Existing Haliotidae records
haliotidae <- RLS_data_year %>%
  filter(family == "Haliotidae") %>%
  distinct(site_code, year, depth)

# Find combinations missing Haliotidae
missing_haliotidae <- all_combos %>%
  anti_join(haliotidae, by = c("site_code", "year", "depth")) %>%
  mutate(
    family = "Haliotidae",
    total = 0
  )
# Add the missing rows
RLS_with_missing_abalone <- bind_rows(RLS_data_year, missing_haliotidae)
#on reviewing this data I feel like all the missing years with abalone data are either a product of a) surveyors being new in 2021, or b)Kieran and claire having these sand sites
#As such, lets not go any further and instead we will keep what we previously had and consider the others as filtered for outliers (we will filter 4 of 7 site_codes anyways due to location or study later other two I think are products of failure to detect / novice surveying)
# CHECK: Review this logic with Isa, I think we should keep 2021 data but not KC data

#Okay, lets come back to year here, and clean the data set. Our goal is to lubridate date, and get density
RLS_data_clean <- RLS_data_year %>%
  mutate(survey_date = ymd(survey_date),
         date_time_survey = ymd_hms(paste(survey_date, hour)), #lubridate date
         site_code = ifelse(site_name == "Swiss Boy", "BMSC24", site_code),
         survey_total = total) %>%
  # just abalone sites (even though we know some have zeros)
  filter(species_name == "Haliotis kamtschatkana") %>% #filter for just the abalone data
  group_by(site_code, year, survey_total)%>%
  select(site_code, site_name, survey_date, date_time_survey, depth, method, species_name, size_class, survey_total) %>% #take relevant variables
  mutate(search = "non_cryptic") #add search column, we will add method later once we dont overright the current method data

# 2026 RLS data from spreadsheet as it is not on the up to date global repository
RLS_new <- read_csv("data-raw/RLS_2026_only_bonus.csv")

#lets check for abalone 0s in the 2026 df
RLS_new_0ab <- RLS_new %>%
  # processing required to get this df into the RLS data format
  filter(Method != 0) %>% # get rid of all method 0's
  slice(2:n()) %>% # cuts the first blank row
  # rename columns
  rename(
    site_code = `Site No.`,
    site_name = `Site Name`, 
    common_name = `Common name`,
    `0` = Inverts,
    species_name = Species,
    method = Method,
    depth = Depth
  )  %>% 
  # Rename columns with spaces
  mutate(
    `15` = na_if(as.character(`15`), "`") #remove back tick in row 
  ) %>%
  mutate(
    across(
      `0`:`400`,
      ~ as.numeric(as.character(.x)) #change all new bonus data to numeric
    )
  )%>%
  mutate(species_name = str_to_sentence(species_name),
         common_name = str_to_sentence(common_name),
         survey_date = dmy(Date),
         date_time_survey = ymd_hms(paste(survey_date, Time))) %>%
  # just abalone
  #filter(species_name == "Haliotis kamtschatkana") %>%
  # Pivot longer for biomass
  pivot_longer(cols = `0`:`400`, names_to = "size_class", values_to = "total") %>% # size class 0 = unsized!!!!
  drop_na(total) %>%
  #filter(total > 0) %>%
  select(-Total)

# All unique site-year-depth combinations
all_combos_2026 <- RLS_new_0ab %>%
  distinct(site_code, depth)

# Existing Haliotidae records
haliotidae_2026 <- RLS_new_0ab %>%
  filter(species_name == "Haliotis kamtschatkana") %>%
  distinct(site_code, depth)

# Find combinations missing Haliotidae
missing_haliotidae_2026 <- all_combos_2026 %>%
  anti_join(haliotidae_2026, by = c("site_code", "depth")) %>%
  mutate(
    species_name = "Haliotis kamtschatkana",
    total = 0
  )
#I'm more confident in these zeros... curious to compare to other years for Hoise south...
#Again, instead of trying to add these rows in using some bind function im just going to mutate two zeros in later in the data set for Hosie south as we know their density will = 0...

#View(RLS_new) #commented out to run all code easily
#starting from the top: we've got a few new issues here a) we got a weird formatting row in row 2 b) we've got M0 data c) the column names are different from the downloaded data... d) the survey blocks aren't added up / together (two surveyers on each side of transect) e) date needs to be lubridated

#lets try to clean this
RLS_new_clean <- RLS_new %>%
  # processing required to get this df into the RLS data format
  filter(Method != 0) %>% # get rid of all method 0's
  slice(2:n()) %>% # cuts the first blank row
  # rename columns
  rename(
    site_code = `Site No.`,
    site_name = `Site Name`, 
    common_name = `Common name`,
    `0` = Inverts,
    species_name = Species,
    method = Method,
    depth = Depth
  )  %>% 
  # Rename columns with spaces
  mutate(species_name = str_to_sentence(species_name),
         common_name = str_to_sentence(common_name),
         survey_date = dmy(Date),
         date_time_survey = ymd_hms(paste(survey_date, Time))) %>%
  # just abalone
  filter(species_name == "Haliotis kamtschatkana") %>%
  mutate(
    `15` = na_if(as.character(`15`), "`") #remove back tick in row 
  ) %>%
  mutate(
    across(
      `0`:`400`,
      ~ as.numeric(as.character(.x)) #change all new bonus data to numeric
    )
  )%>%
  # Pivot longer for biomass
  pivot_longer(cols = `0`:`400`, names_to = "size_class", values_to = "total") %>% # size class 0 = unsized!!!!
  drop_na(total) %>% #removing any n/a from the df
  select(-Total) %>% #dropping previously named "Total"
  # group blocks 1 and 2
  group_by(site_code, site_name, survey_date, date_time_survey, depth, method, species_name, common_name, size_class) %>%
  summarise(survey_total = sum(total)) %>% # sum blocks 1 and 2 (these are either side of the transect)
  ungroup() %>%
  mutate(size_class = as.numeric(size_class),
         # correct for rectangle area
         survey_den = case_when(method == 1 ~ survey_total/500,
                                method == 2 ~ survey_total/100)) %>% #dont end up using this
  as.data.frame() %>%
  select(site_code, site_name, survey_date, date_time_survey, depth, method, species_name, size_class, survey_total) %>%
  mutate(year = year(survey_date)) %>% #add year column
  mutate(search = "non_cryptic") #add search column, again will include method later once no need for pre-exisiting method data referencing fish or invert counts in the df



#lets merge the data sets using rbind
RLS_full <- rbind(RLS_new_clean, RLS_data_clean)

#View(RLS_full) #commented out to run all code easily
#notice survey density isn't on the data frame, the manipulation above only worked per size class. Need to think of a way to capture for all size classes across a survey site and depth.
#If we don't care about depth we will then sum those together and divide by the appropriate area

RLS_data_totals <- RLS_full %>%
  group_by(site_name, site_code, year, depth)%>%
  summarise(total = sum(survey_total)) #lets take just what we want and get a sum total across all size classes of abalone

#now build code that checks if year and site_code occur more than once, if that is the case sum the total and divide by 200, otherwise just divide the total by 100. Reason for this is some sites are surveyed across two depths at the same time in the same year. Therefore, we need to sum their counts and divide by twice the area surveyed in a single survey to get a site level mean density
RLS_data_means <- RLS_data_totals %>%
  group_by(site_name, site_code, year) %>%
  summarise(
    site_mean_density = if (n() == 2) {
      sum(total, na.rm = TRUE) / 200
    } else {
      first(total) / 100
    },
    .groups = "drop"
  )

#TO DO check against how many should drop from RLS_data_totals (have an expectation based on knowledge of sites and years for a CYU)

#Now we create rows that adds zero for BMSC26, year = 26, density =0
#remove sites that aren't repeated, Kieran cox studies on subtidal kelp beds, or broken group sites
unique(RLS_data_means$site_code) #see the ones I need to remove

RLS_data_means_filtered <- RLS_data_means %>%
  rename(site = "site_name") %>%
  add_row(
    site = "Hosie South",
    site_code = "BMSC26",
    year = 2026,
    site_mean_density = 0)%>%
  add_row(site = "Hosie South",
          site_code = "BMSC26",
          year = 2026,
          site_mean_density = 0)%>% #we just added the relevant 0s found in missing_haliotidae_2026
  filter(!site_code %in% c(
    "BMKC2",
    "BMSC31",
    "BMSC32",
    "BMSC33",
    "KCCA11",
    "KCCA13",
    "KCCA19" #, #all this data is from a separate study in 2023 that had differing objectives than just rocky reef surveys
   # "BMSC13",
   # "BMSC14",
   # "BMSC15",
   # "BMSC16",
   # "BMSC17",
   # "BMSC18" #all this data is sites in the broken group but lets keep this for now
  ))%>%
  mutate(method = "par_belt")%>%
  mutate(search = "non_cryptic")%>%
  mutate(site = case_when(str_detect(site, "Aguilar Point") ~ "aguilar_point", #changing names relevant to outplanting information.
                          str_detect(site, "Eagle Bay") ~ "scotts_bay",
                          str_detect(site, "Goby Town") ~ "grappler",
                          TRUE ~ site)) 
#taylor or dodgers could be included as the ed_king_se most adjacent site here but it doesn't quite feel correct

#okay lets check if our removing sites worked
unique(RLS_data_means_filtered$site_code)
#Okay yay! We got there, lets save that as processed data
#readr::write_csv(RLS_data_means_filtered, file = "data-processed/RLS_data_2021_2026.csv") #to save the processed data (commented out so that when running the code you dont save something accidentally that you may have changed)

##lets get an average for each site across all years 2021-2026
RLS_average_density <- RLS_data_means_filtered %>%
  group_by(site_code) %>%
  summarise(overall_mean_density = mean(site_mean_density), #calculate density average across all years of study
            number_years = n_distinct(year), #tell us how many years go into the average calculation
            .groups = "drop")
#hmm for some reason many sites have less years sampled than expected. Look into this tomorrow