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
  mutate(year = if_else(is.na(year), 1994, year)) #add year in places when they weren't directly pasted into the raw file but they were in the document titled 1994

  
#we are going to leave the size information for JW right now as really dont have much of that...
#no depth for 1994 - perhaps it was assumed to be the same as previous survey?
#also note depth makes major chage throughout program

View(JW_data_clean)
#readr::write_csv(JW_data_clean, file = "data-processed/JW_data_1988_2022.csv") #to save the clean data (commented out so that when running the code you dont save something accidentally that you may have changed)

#TO DO CONFIRM METHODS

#Christine Hansen data -------------------------------------------------------------------
#in the future I hope to have this raw data...
#for now, to get this data we have to use meta digitize 
#CH_data_shallow <- metaDigitise("/Users/carterburtlake/Documents/RProjects/MScChapter2/figures")
#load in shallow data

#CH_data_deep <- metaDigitise("/Users/carterburtlake/Documents/RProjects/MScChapter2/figures")
#load in deep data

#readr::write_csv(CH_data_deep, file = "data-raw/CH_data_2008.csv") #to save the raw data

CH_data_deep <- read.csv("data-raw/CH_data_2008.csv")
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
                          str_detect(group_id, "helby") ~ "helby_sw",
                          TRUE ~ "scotts_bay")) %>% #give them independent site names outside of depth and reproductive status
  add_row(study = "ChristineHansen", mean = 0, se = 0, year = "2008", depth_m = "7.5", reproductive_status = "reproductive", site = "helby_sw")%>%
  add_row(study = "ChristineHansen", mean = 0, se = 0, year = "2008", depth_m = "2.5", reproductive_status = "reproductive", site = "helby_sw")%>%
  add_row(study = "ChristineHansen", mean = 0, se = 0, year = "2008", depth_m = "7.5", reproductive_status = "reproductive", site = "grappler")%>%
  add_row(study = "ChristineHansen", mean = 0, se = 0, year = "2008", depth_m = "7.5", reproductive_status = "immature", site = "grappler")%>%
  add_row(study = "ChristineHansen", mean = 0, se = 0, year = "2008", depth_m = "2.5", reproductive_status = "reproductive", site = "grappler")%>%
  add_row(study = "ChristineHansen", mean = 0, se = 0, year = "2008", depth_m = "2.5", reproductive_status = "immature", site = "grappler")%>% #had to add all the goby town data where they searched for abalone but found none during their surveys. Same is also true from "reproductive" abalone
  mutate(method = "par_perp_quad")%>% #method used is similar to breen
  mutate(search = "cryptic") #every fourth quadrat they flipped rocks
#note that year, and depth are characters... we will probably need them as numeric and dates at some point so lets change that now
is.character(CH_data_clean$year) #TRUE
is.character(CH_data_clean$depth_m) #TRUE


CH_data_mean <- CH_data_clean%>%
  group_by(site, year, method, search)%>%
  summarise(site_mean_density = sum(mean, na.rm = TRUE))
#add the means across the whole site regardless of depth or size class / reproductive status 
#this is processed data and should be saved as such

#readr::write_csv(CH_data_mean, file = "data-processed/CH_data_2008.csv") #to save the processed data (commented out so that when running the code you dont save something accidentally that you may have changed)

###FUTURE: NEED TO GET DEPTH DATA IN THIS DF
#DFO kelp data ------------------------------------------------------------------------
DFO_kelp_data <- read.csv("data-raw/2021-2025_DFO_Kelp_Abalone.csv")

#hmm this data set looks different... I need to get year and site into this data frame...
#From the methods I know the area covered on each survey, note that some comments or notes might describe missing transect lines or quadrats
#for the purposes of this project I'm going to leave this data set out

#Storm Coast data --------------------------------------------------------------------------------
SC_data <- read.csv("data-raw/StormCoast_Abalone_2025.csv")

View(SC_data)
is.character(SC_data$date)
###review dixon inside on datasheets -> for now filter the site out
#TO DO: fix column names, get dates to year, go from sizes to density per site
SC_data <- SC_data %>%
  clean_names() #make column names easier to work with
SC_data_clean <- SC_data%>%
  filter(quadrat != "NA") %>% #remove a bunch of blank data, now this matches our raw df...
  filter(site_id != "Dixon Inside") %>% #remove site where we are uncertain of dead or alive counts - REVIEW PAPER DATA TO VALIDATE THIS 
  mutate(date = ymd(date), year = year(date)) %>% #change date to a date using lubridate and add a year column
  mutate(method = "perp_quad") %>% #method used 
  mutate(search = "non_cryptic") %>% #did not flip over rocks in search
  mutate(site_id = case_when(str_detect(transect_id, "EKE02") ~ "ed_king_se", #changing names relevant to outplanting information.
                          str_detect(transect_id, "EGB01") ~ "scotts_bay",
                          str_detect(transect_id, "ELI01") ~ "ellis",
                          str_detect(transect_id, "GOT01") ~ "grappler",
                          str_detect(transect_id, "HSW01") ~ "helby_sw",
                          str_detect(transect_id, "SAN02") ~ "sandford_sw", #this site might be less accurate -> to review in CH thesis
                          str_detect(transect_id, "AGU01") ~ "aguilar_point",
                          TRUE ~ site_id)) 
  
 

unique(SC_data_clean$alive_or_dead) #make sure next alive filter will correctly get all data
#Nevermind, this wont work as it will impact our number of quadrat count...
#instead lets find a total alive value per site and total quadrat per site and join them
SC_data_alive <- SC_data_clean %>%
  group_by(site_id, year, method, search)%>%
  summarise(total_alive = sum(alive_or_dead == "alive")) #first we get number of alive per site

SC_data_quad <- SC_data_clean %>%
  group_by(site_id, year, method, search)%>%
  summarise(n_quadrats = n_distinct(quadrat)) #then we get number of quads per site

SC_data_join <- left_join(SC_data_alive, SC_data_quad, by = c("site_id", "year", "method", "search")) %>%
  mutate(site_mean_density = total_alive/n_quadrats) %>% #finally we join them together and calculate the site mean using the number of alive individuals at the site by the number of 1m2 surveyed quads
  rename(site = site_id) #rename to match other data sets

#notice we surveyed 56 sites in this effort but we filtered out dixon inside so 55 observations provides a good moment to CYU
  
#lets save this (commented out so that when running the code you dont save something accidentally that you may have changed)
#readr::write_csv(SC_data_join, file = "data-processed/SC_data_2025.csv") #to save the processed data (commented out so that when running the code you dont save something accidentally that you may have changed)
#RLS data ----------------------------------------------------------------------------------------
RLS_data <- read.csv("data-raw/RLS_2025.csv") 

#lets review the data
View(RLS_data)

#lets try to filter for sites over years that dont have abalone 
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
#As such, lets not go any further and instead we will keep what we previously had and consider the others as filtered for outliers (we will filter 5 of 7 site_codes anyways due to location or study later other two I think are products of failure to detect / novice surveying)

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
RLS_new <- read_csv("data-raw/RLS_2026_only.csv")

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

View(RLS_new)
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

View(RLS_full)
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
    "KCCA19", #all this data is from a separate study in 2023 that had differing objectives than just rocky reef surveys
    "BMSC13",
    "BMSC14",
    "BMSC15",
    "BMSC16",
    "BMSC17",
    "BMSC18" #all this data is sites in the broken group 
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

#Join abalone density-----------------------------------------------------------------------------
#okay now I need to join the data into one df. Lets review them all
#get them all down to 5 variables as listed below
View(JW_data_clean)
JW_join <- JW_data_clean%>%
  select(c("site", "year", "site_mean_density", "method", "search"))
View(CH_data_mean) #this one is good
View(SC_data_join)
SC_join <- SC_data_join%>%
  select(c("site", "year", "site_mean_density", "method", "search"))
View(RLS_data_means_filtered)
RLS_join <- RLS_data_means_filtered%>%
  select(c("site", "year", "site_mean_density", "method", "search")) #remove site code

#I want the variables site, year, mean_site_density, method, and search across all studies in one df
#before combining I first need to make a decision on site overlap (which sites are considered to be the same survey location or not)
#Struggling on decision criteria... lets try two versions?

#lets just bind them with all these different site names and figure out that detail later
abalone_density <- rbind(JW_join, CH_data_mean, SC_join, RLS_join)
#CYU: previously had 84+4+110+55 observations 
84+4+110+55 #matches observations of bound df at 253

View(abalone_density) #I need to add outplant details to this df 

#Outplant information -----------------------------------------------------------------------------
outplant <- read.csv("data-raw/outplant_read_raw.csv")

View(outplant)
#wow this is ugly. For the purposes of this project I'm just going to take the outplant information, not experimentation that resulted in outplants as verifying location of effective outplant during these was not the primary goal 

clean_outplant <- clean_names(outplant) %>% #lets clean this up, start with column names
  filter(str_starts(objective, "BHCAP")) %>% #just take the outplants through BHCAP
  mutate(year = as.integer(str_extract(date, "\\d{4}"))) %>% #lets get the year out of this unorganized weird date column
  mutate(number_outplanted = as.numeric(str_remove_all(amount, ","))) %>% # lets get rid of the commas in our numeric data for number of outplants
  #last thing to do here is to change the outplant names to the specific location of outplant that can correspond to our joined df about abalone density
  mutate(outplant_site = case_when(str_detect(outplant_site, "Scotts Bay") ~ "scotts_bay", 
                          str_detect(outplant_site, "Helby Island") ~ "helby_sw",
                          str_detect(outplant_site, "Grappler Inlet") ~ "grappler",
                          str_detect(outplant_site, "Edward King Island") ~ "ed_king_se",
                          str_detect(outplant_site, "Aguilar Point") ~ "aguilar_point",
                          str_detect(outplant_site, "Sandford Island") ~ "sandford_sw")) 
  

#Join all data -------------------------------------------------------------------------------------
  


#brief visualization
ggplot(data = JW_join, aes(x = year, y = site_mean_density, colour = site)) +
  geom_point()
