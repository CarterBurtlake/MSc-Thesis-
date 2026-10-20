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
  mutate(Density = total_abalone / Area) %>% #calc density
  mutate(method = "par_quad",  #add columns to highlight method, source, and search type
         source = "DFO",
         search = "non_cryptic")

#cant run this twice, must clean the environment if so to make changes
DFO_abalone_density <- clean_names(DFO_abalone_density) %>% #tidy the column names
  rename(site = site_name) %>% #rename site to match all dfs
  rename(site_mean_density = density) %>% #rename to match all dfs
  select(site,year,site_mean_density,source,search,method) %>% #remove columns to match all dfs
  mutate(site = str_to_lower(str_replace_all(site, " ", "_")), #tidy the site names
         site = str_replace_all(site, "'", "")) 

#lets get an average for each site across all years 2021-2025
DFO_average_density <- DFO_abalone_density %>%
  filter(!site %in% c("taylor_island_2021", "flemming_112_2021", "flemming_114")) %>% #remove sites that were only done in 2021
  # CHECK: Confirm we should remove or keep these? I think we should remove in this specific case as the sites were seen as unsuitable
  group_by(site) %>%
  summarise(average_density = mean(site_mean_density), #calculate density average across all years of study
            number_years = n_distinct(year), #tell us how many years go into the average calculation
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

RLS_data_totals_missing <- RLS_data_totals %>%
  ungroup()%>%
  add_row(site_name = "Dixon Inside", #add survey data that says none were found here in 2021 (removed when filtered for abalone)
          site_code = "BMSC22",
          year = 2021,
          total = 0) %>%
  add_row(site_name = "Raymond Kelp Rock", #add survey data that says none were found here in 2021 (removed when filtered for abalone)
          site_code = "BMSC15",
          year = 2021,
          total = 0) %>%
  add_row(site_name = "Hosie South", #add survey data that says none were found here in 2026 (removed when filtered for abalone)
          site_code = "BMSC26",
          year = 2026,
          total = 0) %>%
  add_row(site_name = "Wizard Island South", #add survey data that says none were found at this depth in 2022 (removed when filtered for abalone, not in 2022 one survey found abalone but there are two surveys so if I want the calculation to appropriately score density below I need to add another row that has a 0 total)
          site_code = "BMSC11",
          year = 2022,
          total = 0)
#now build code that checks if year and site_code occur more than once, if that is the case sum the total and divide by 200, otherwise just divide the total by 100. Reason for this is some sites are surveyed across two depths at the same time in the same year. Therefore, we need to sum their counts and divide by twice the area surveyed in a single survey to get a site level mean density
RLS_data_means <- RLS_data_totals_missing %>%
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

#tidy names to match combined df of all average abalone density
RLS_abalone_density <- clean_names(RLS_data_means_filtered) %>% #tidy the column names
  select(site, year, site_mean_density, method, search) %>% #remove site code to match other dfs
  mutate(site = str_to_lower(str_replace_all(site, " ", "_")), #tidy the site names
         site = str_replace_all(site, "'", ""),
         source = "RLS") 

#okay lets check if our removing sites worked
unique(RLS_data_means_filtered$site_code)
#Okay yay! We got there, lets save that as processed data
#readr::write_csv(RLS_data_means_filtered, file = "data-processed/RLS_data_2021_2026.csv") #to save the processed data (commented out so that when running the code you dont save something accidentally that you may have changed)

##lets get an average for each site across all years 2021-2026
RLS_average_density <- RLS_data_means_filtered %>%
  group_by(site_code) %>%
  summarise(average_density = mean(site_mean_density), #calculate density average across all years of study
            number_years = n_distinct(year), #tell us how many years go into the average calculation
            .groups = "drop")
#Reviewed all sites and all have the correct number of years of data contributing to their average
# CHECK: worth investigating density by depth 
# CHECK: see if its worth adding KCCA data...

# Parks data ----
#load 2007-2012 data
Parks_data_2007 <- read.csv("data-raw/PRNPR-Northern-Abalone-2007-2012.csv") %>% #size data
  rename(site = `Location`)
Parks_data_2007_sites <- read.csv("data-raw/PRNPR-Northern-Abalone-2007-2012_sitedata.csv") #site data
  

#clean size dataset and obtain total_abalone per site
Parks_clean_data_2007 <- clean_names(Parks_data_2007) %>%
  filter(!comment == "dead") %>%
  group_by(site, year) %>%
  summarise(total_abalone = n()) #add up the number of measured abalone for a year and site

#clean site dataset and obtain site level survey_area
Parks_clean_sites_2007 <- clean_names(Parks_data_2007_sites) %>%
  group_by(site, year)%>%
  summarise(survey_total_transect = (n_distinct(transect)), #number of transects at each site in a given year
            survey_area = survey_total_transect * 3.75) #number of transects * area surveyed in a single transect

#join datasets together and get density
Parks_abalone_2007 <- Parks_clean_sites_2007 %>%
  left_join(Parks_clean_data_2007 %>%
              select(site, year, total_abalone), #add the total ab data to the site data
            by = c("site", "year")) %>% #combining the dataset by year and site
  mutate(total_abalone = replace_na(total_abalone, 0)) %>% #N/A -> 0s for years when abalone wasnt seen at a site
  mutate(density = total_abalone / survey_area) #calc density

#average density
average_Parks_abalone_2007 <- Parks_abalone_2007 %>%
  group_by(site) %>%
  summarise(mean_density = mean(density), #calculate density average across all years of study
            number_years = n_distinct(year), #tell us how many years go into the average calculation
            .groups = "drop")

#need to join with the 2013-2016 data and KRI data to get the complete average density across these sites

#load 2013-2016 data
Parks_data_2013 <- read.csv("data-raw/PRNPR-Northern_Abalone-2013-2016.csv") #size data
Parks_data_2013_sites <- read.csv("data-raw/PRNPR-Northern_Abalone-2013-2016_sitedata.csv") #site data

#clean size dataset and obtain total_abalone per site
Parks_clean_data_2013 <- clean_names(Parks_data_2013) %>%
  group_by(site, year) %>%
  summarise(total_abalone = n()) #add up the number of distinct measured abalone for a year and site

#clean site dataset and obtain site level survey_area
Parks_clean_sites_2013 <- clean_names(Parks_data_2013_sites) %>%
  mutate(site = trimws(site)) %>% #one of the 2016 transects had a space after the name "Dodd " so use trimws to fix
  group_by(site, year)%>%
  summarise(survey_total_transect = (n_distinct(transect)), #number of transects at each site in a given year
            survey_area = survey_total_transect * 60) #number of transects * area surveyed in a single transect

#join datasets together and get density
Parks_abalone_2013 <- Parks_clean_sites_2013 %>%
  left_join(Parks_clean_data_2013 %>%
              select(site, year, total_abalone), #add the total ab data to the site data
            by = c("site", "year")) %>% #combining the dataset by year and site
  mutate(density = total_abalone / survey_area) #calc density

#average density
average_Parks_abalone_2013 <- Parks_abalone_2013 %>%
  group_by(site) %>%
  summarise(mean_density = mean(density), #calculate density average across all years of study
            number_years = n_distinct(year), #tell us how many years go into the average calculation
            .groups = "drop")

#now we need to merge all parks density data to calc a total parks average density data (but we dont really need this for thesis)
#lets merge the data sets using rbind
Parks_full_abalone <- rbind(Parks_abalone_2013, Parks_abalone_2007) %>%
  mutate(site = trimws(site)) %>% #effingham name has white space, remove so that they are seen as same site
  mutate(site = case_when(str_detect(site, "Turret") ~ "Turret Island", #rename turret to match 2007-2012 data
                          TRUE ~ site))

#tidy names to match combined df of all average abalone density
Parks_abalone_density <- clean_names(Parks_full_abalone) %>% #tidy the column names
  rename(site_mean_density = density) %>% #rename to match other dfs naming convention
  select(site, year, site_mean_density) %>% #remove site code to match other dfs
  mutate(site = str_to_lower(str_replace_all(site, " ", "_")), #tidy the site names
         site = str_replace_all(site, "'", ""),
         source = "Parks",
         method = "par_belt",
         search = "non_cryptic") 

#average density
average_Parks_abalone <- Parks_full_abalone %>%
  group_by(site) %>%
  summarise(mean_density = mean(density), #calculate density average across all years of study
            number_years = n_distinct(year), #tell us how many years go into the average calculation
            .groups = "drop")
# KRI data ----
#load data
KRI_data<- read.csv("data-raw/BaumLab_BGI_abalonedata_2021-2025.csv") 

#clean data set
KRI_data_clean <- clean_names(KRI_data) 

#clean dates
KRI_data_clean <- KRI_data_clean %>%
  mutate(date = as.Date(date, format = "%d.%m.%Y"))%>%
  mutate(year = year(date)) %>% #get year column
  mutate(site = case_when(str_detect(site, "Turret") ~ "Turret Island", #rename to ensure all names in df match
                          str_detect(site, "Effingham") ~ "Effingham Bay", #rename to match Parks data
                          str_detect(site, "Hand") ~ "Hand Island", #rename to match Parks data
                          str_detect(site, "Nettle") ~ "Nettle", #rename to match Parks data
                          str_detect(site, "Harb") ~ "Island Harbour", #rename to match Parks data
                          TRUE ~ site)) 

#get survey area 
KRI_data_survey_area <- KRI_data_clean %>%
  group_by(site, year)%>%
  summarise(survey_total_transect = (n_distinct(transect)), #number of transects at each site in a given year
            survey_area = survey_total_transect * 60) #number of transects * area surveyed in a single transect

KRI_abalone_total <- KRI_data_clean %>%
  filter(size_cm > 0)%>% #remove rows where size = 0 as there were no abalone
  group_by(site, year) %>%
  summarise(total_abalone = n()) #total the number of abalone by totaling the number of rows that contain both year and site info CYU if just look at the df from the filter we see there are 68 observations, then the total abalone count across all sites and years add up to that so we know each row that = 68 so each size is tallied is accounted for -> total abs
  
KRI_abalone_totals_with_0s <- KRI_abalone_total%>%  #add 0s for sites that found no abalone
  ungroup() %>%  # Remove grouping before adding rows
  add_row(site = "Effingham Bay", 
          year = 2021,
          total_abalone = 0) %>%
  add_row(site = "Dodd", 
          year = 2021,
          total_abalone = 0) %>%
  add_row(site = "Hand Island", 
          year = 2025,
          total_abalone = 0) %>%
  add_row(site = "Island Harbour", 
          year = 2021,
          total_abalone = 0) %>%
  add_row(site = "Island Harbour", 
          year = 2025,
          total_abalone = 0)
  
KRI_density <- KRI_data_survey_area %>%
  left_join(KRI_abalone_totals_with_0s %>%
              select(site, year, total_abalone), #add the total ab data to the site data
            by = c("site", "year")) %>% #combining the dataset by year and site
  mutate(density = total_abalone / survey_area) #calc density

#tidy names to match combined df of all average abalone density
KRI_abalone_density <- clean_names(KRI_density) %>% #tidy the column names
  rename(site_mean_density = density) %>% #rename to match other dfs naming convention
  select(site, year, site_mean_density) %>% #remove site code to match other dfs
  mutate(site = str_to_lower(str_replace_all(site, " ", "_")), #tidy the site names
         site = str_replace_all(site, "'", ""),
         source = "KRI",
         method = "par_belt",
         search = "non_cryptic") 

#average density
average_KRI_abalone <- KRI_density %>%
  group_by(site) %>%
  summarise(mean_density = mean(density), #calculate density average across all years of study
            number_years = n_distinct(year), #tell us how many years go into the average calculation
            .groups = "drop")

# Parks + KRI data ----
#lets merge all the density data for the Parks sites
BGI_abalone_density <- rbind(Parks_abalone_density, KRI_abalone_density)

#average density across all the years (just for fun)
average_BrokenGroup_data <- BGI_abalone_density %>%
  group_by(site) %>%
  summarise(average_density = mean(site_mean_density), #calculate density average across all years of study
            number_years = n_distinct(year), #tell us how many years go into the average calculation
            .groups = "drop")

#plots sites through the years 
test <- ggplot(data = BGI_abalone_density, aes(x = as.numeric(year), y = site_mean_density, colour = site)) +
  geom_point()+
  geom_line()+
  labs(x = "Year", y = "Abalone density (" ~m^-2*")", colour = "Site") + #add titles
  theme_classic()

test
#notice that the 2013-2016 surveys had the largest search effort and seem to have the most variability in the data. 


# Storm Coast data ----
#load data
SC_data <- read.csv("data-raw/StormCoast_Abalone_2025.csv")

View(SC_data)#commented out to run all code easily
is.character(SC_data$date)
# CHECK: review dixon inside on datasheets -> for now filter the site out

#clean data
SC_data <- SC_data %>%
  clean_names() #make column names easier to work with

SC_data_clean <- SC_data%>%
  filter(quadrat != "NA") %>% #remove a bunch of blank data, now this matches our raw df... doesnt remove na's in rows without abalone though which is perfect
  filter(site_id != "Dixon Inside") %>% #remove site where we are uncertain of dead or alive counts - REVIEW PAPER DATA TO VALIDATE THIS 
  filter(site_id != "Flemming 114") %>% #remove site where we are uncertain of dead or alive counts - REVIEW PAPER DATA TO VALIDATE THIS 
  mutate(date = ymd(date), year = year(date)) %>% #change date to a date using lubridate and add a year column
  mutate(method = "perp_quad") %>% #method used 
  mutate(search = "non_cryptic") %>% #did not flip over rocks in search
  mutate(source = "Storm_Coast") %>% #add to match df convention
  mutate(site_id = case_when(str_detect(transect_id, "EKE02") ~ "ed_king_se", #changing names relevant to outplanting information.
                             str_detect(transect_id, "EGB01") ~ "scotts_bay",
                             str_detect(transect_id, "ELI01") ~ "ellis",
                             str_detect(transect_id, "GOT01") ~ "grappler",
                             str_detect(transect_id, "HSW01") ~ "helby_sw",
                             str_detect(transect_id, "SAN02") ~ "sandford_sw", #this site might be less accurate -> to review in CH thesis
                             str_detect(transect_id, "AGU01") ~ "aguilar_point",
                             TRUE ~ site_id)) 




#find a total alive value per site and total quadrat per site and join them
SC_data_alive <- SC_data_clean %>%
  group_by(site_id, year, method, search,source)%>%
  summarise(total_alive = sum(alive_or_dead == "alive")) #first we get number of alive per site

SC_data_quad <- SC_data_clean %>%
  group_by(site_id, year, method, search, source)%>%
  summarise(n_quadrats = n_distinct(quadrat)) #then we get number of quads per site

SC_data_join <- left_join(SC_data_alive, SC_data_quad, by = c("site_id", "year", "method", "search", "source")) %>%
  mutate(site_mean_density = total_alive/n_quadrats) %>% #finally we join them together and calculate the site mean using the number of alive individuals at the site by the number of 1m2 surveyed quads
  rename(site = site_id) #rename to match other data sets

unique((SC_data_join$site))
#notice we surveyed 57 sites in this effort but we filtered out dixon inside and flemming 114 so 55 observations provides a good moment to CYU
#CHECK : get paper data and validate

#tidy names to match combined df of all average abalone density
SC_abalone_density <- clean_names(SC_data_join) %>% #tidy the column names
  select(site, year, site_mean_density,method,source,search) %>% #remove other variables to match dfs 
  mutate(site = str_to_lower(str_replace_all(site, " ", "_")), #tidy the site names
         site = str_replace_all(site, "'", ""))

#lets save this (commented out so that when running the code you dont save something accidentally that you may have changed)
#readr::write_csv(SC_data_join, file = "data-processed/SC_data_2025.csv") #to save the processed data (commented out so that when running the code you dont save something accidentally that you may have changed)

# Jane Watson data ----
#load data  
JW_data <- read.csv("data-raw/AbaloneJaneWatsonCombined89_22.csv")
#all sites accounted for (84)

#clean data set
JW_data_clean <- clean_names(JW_data) %>%
  mutate( #tidy site names
    code = as.integer(str_extract(site, "^\\d+")),
    site = case_when(
      code == 1 ~ "ship_islands",
      code == 2 ~ "ed_king_sw", #call relevant to outplant location
      code == 3 ~ "ed_king_nw", #call relevant to outplant location
      code == 4 ~ "seppings",
      code == 5 ~ "cape_beale",
      code == 6 ~ "lawton_point",
      code == 7 ~ "whittlestone",
      code == 8 ~ "execution_rock",
      code == 9 ~ "self_point", #this sites gps location is near effingham inlet?
      code == 10 ~ "helby_ne",  #call relevant to outplant location
      code == 11 ~ "scotts_bay",
      code == 12 ~ "aguilar_point",
      code == 13 ~ "cia_rock",
      code == 14 ~ "ed_king_se", #for the purposes of this project I am renaming this site to align with the outplant location name (it is the closest site in this study)
      code == 15 ~ "kirby_point",
      code == 16 ~ "village_bay",
      code == 17 ~ "blowhole",
      code == 18 ~ "prasiola",
      code == 19 ~ "wizard",
      code == 20 ~ "grappler",
      code == 21 ~ "ed_king_?", #seems like the same place as ed_king_nw
      TRUE ~ site
    )
  ) %>%
  select(-code)%>%
  mutate(density_sample_size = str_remove(density_sample_size, "^n="))%>% #remove n= from rows in density sample size
  mutate(site_mean_density = mean/0.49)%>% #okay so quadrats are 0.7m x 0.7m = 0.49m^2, so, to get abalone mean in quadrats across 1m2 we must do mean abalone in 0.49m2 / 0.49m2 = abalone per m2
  mutate(depth_m = depth_ft * 0.3048)%>%
  mutate(date = mdy(date), year = year(date))%>% #changes all the date formats using lubridate so that the second call of "year" can be applied to make a new column with just year
  mutate(method = "par_quad") %>% #add the method of search used (in this case it was parallel to shore at a depth between 5-12m with quads a certain number of fin kicks away (Watson & Estes, 2011))
  mutate(search = "non_cryptic") %>% #add the search type
  mutate(source = "Jane_Watson") %>% #add the source 
  mutate(year = if_else(is.na(year), 1994, year)) #add year in places when they weren't directly pasted into the raw file but they were in the document titled 1994

# CHECK: confirm with Jane that the methods used for these surveys different than the wizard and taylor dataset?
# CHECK: important to stress that I am just taking Janes mean data and we could recalculate these data
# CHECK: there might be additional abalone in this dataset that I haven't accessed wherein they used independent quadrats to sample abalone and urchin sizes
#we are going to leave the size information for JW right now as really dont have much of that...
#no depth for 1994 - perhaps it was assumed to be the same as previous survey?
#also note depth makes major change throughout program

#tidy names to match combined df of all average abalone density
JW_abalone_density_7year <- clean_names(JW_data_clean) %>% #tidy the column names
  select(site, year, site_mean_density,method,source,search) %>% #remove other variables to match dfs 
  mutate(site = str_to_lower(str_replace_all(site, " ", "_")), #tidy the site names
         site = str_replace_all(site, "'", ""))

#View(JW_data_clean) #commented out to run all code easily
#readr::write_csv(JW_data_clean, file = "data-processed/JW_data_1988_2022.csv") #to save the clean data (commented out so that when running the code you dont save something accidentally that you may have changed)

#lets get an average for each site across all years 1989-2022
JW_average_density <- JW_abalone_density_7year %>%
  filter(!site == "ed_king_?") %>% #remove site with no data?
  # CHECK: what is this data???
  group_by(site) %>%
  summarise(average_density = mean(site_mean_density), #calculate density average across all years of study
            number_years = n_distinct(year), #tell us how many years go into the average calculation
            .groups = "drop")

#load annual data set
JW_taylor_data <- read.csv("data-raw/JW_Taylor_Annual.csv")
JW_wizard_data <- read.csv("data-raw/JW_Wizard_Annual.csv")

#clean taylor
JW_taylor_data_clean <- clean_names(JW_taylor_data)%>%
  select(year,abalone) %>% #take column that totals abalone and the year it happened
  filter(str_detect(year, regex("MEAN", ignore_case = TRUE))) %>% #could do this for variance and SD
  mutate(year = as.numeric(str_extract(year, "\\d{2}$")), #extract number at end of name, then add appropriate year
         year = if_else(year <= 26, year + 2000, year + 1900))%>%
  rename(mean_abalone = abalone) %>% #rename to be clear
  mutate(survey_area = 20, #each survey has 5 transect lines of 10m that is surveyed 1m on either side ((2x10)x5) = 100m of survey area *but since pulling means already just need to divide by indivudal transect and get the same value (i.e, total abalone counted across 5 transects / 100m2 = mean abalone counted across all 5 surveys / 20m2)
         site = "taylor") 

#density taylor
JW_taylor_density <- JW_taylor_data_clean%>%
  filter(!year == "1991")%>% #remove year that has no data
  mutate(mean_abalone = as.numeric(mean_abalone))%>% #make these mean values numeric so can calculate density
  mutate(site_mean_density = mean_abalone/survey_area) #calc density per year 
  
#do the same with Wizard
JW_wizard_data_clean <- clean_names(JW_wizard_data)%>%
  select(year,abalone) %>% #take column that totals abalone and the year it happened
  filter(str_detect(year, regex("MEAN", ignore_case = TRUE))) %>% #could do this for variance and SD
  mutate(year = as.numeric(str_extract(year, "\\d{2}$")), #extract number at end of name, then add appropriate year
         year = if_else(year <= 26, year + 2000, year + 1900))%>%
  rename(mean_abalone = abalone) %>% #rename to be clear
  mutate(survey_area = 20,
         site = "wizard") 

#density wizard
JW_wizard_density <- JW_wizard_data_clean%>%
  mutate(site_mean_density = mean_abalone/survey_area) #calc density per year 

#annual surveys combined and get in same format as other dfs
JW_annual_density <- rbind(JW_wizard_density, JW_taylor_density) %>%
  select(site,year,site_mean_density) %>% #remove columns to keep dfs consistent
  mutate(search = "non_cryptic",
         source = "Jane_Watson",
         method = "perp_belt")

#combine all Jane Watson data
JW_abalone_density <- rbind(JW_abalone_density_7year, JW_annual_density)

# SubSci data ----

# Combined sources ----
abalone_density <- rbind(JW_abalone_density, DFO_abalone_density, RLS_abalone_density, KRI_abalone_density, SC_abalone_density)

#check site alignment
unique(abalone_density$site)
#figure out 1 year DFO sites (ask sandie if decided to be unsuitable habitat)
#figure out naming convention in map for 100m proximity expcept if aspect changes (have to maybe filter wizard, clutus, folger, flemming, taylor/dodger/edking se complex, Sanford, ed_king all)

#could do all this with pychno data across all data sets except subsci 2024