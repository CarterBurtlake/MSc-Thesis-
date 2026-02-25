# For Carter

# load packages
library(tidyverse)


# RLS data from the website!
test <- read.csv("Data/mobile_macroinvertebrate_abundance.csv")
rls <- read_csv("Data/mobile_macroinvertebrate_abundance.csv",
                 show_col_types = FALSE) %>%
  as.data.frame() %>%
  mutate(survey_date = ymd(survey_date),
         date_time_survey = ymd_hms(paste(survey_date, hour)),
         site_code = ifelse(site_name == "Swiss Boy", "BMSC24", site_code),
         # correct for rectangle area
         survey_den = case_when(method == 1 ~ total/500,
                                method == 2 ~ total/100),
         survey_total = total) %>%
  # just abalone
  filter(species_name == "Haliotis kamtschatkana") %>%
 # filter(month(survey_date) == 4 | month(survey_date) == 5) %>% # can use to filter for just the blitz data
  select(site_code, site_name, survey_date, date_time_survey, depth, method, species_name, size_class, survey_total)

# RLS data from spreadsheet
kelp_rls <- read_csv("Data/RLS_KCCA_2022.csv") %>%
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
  drop_na(total) %>%
  filter(total > 0) %>%
  select(-Total) %>%
  # group blocks 1 and 2
  group_by(site_code, site_name, survey_date, date_time_survey, depth, method, species_name, common_name, size_class) %>%
  summarise(survey_total = sum(total)) %>% # sum blocks 1 and 2
  ungroup() %>%
  mutate(size_class = as.numeric(size_class),
         # correct for rectangle area
         survey_den = case_when(method == 1 ~ survey_total/500,
                                method == 2 ~ survey_total/100)) %>%
  as.data.frame() %>%
  select(site_code, site_name, survey_date, date_time_survey, depth, method, species_name, size_class, survey_total)

# 2025 RLS data from spreadsheet
new_rls <- read_csv("Data/RLS_2025_Data.csv") %>%
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
  drop_na(total) %>%
  filter(total > 0) %>%
  select(-Total) %>%
  # group blocks 1 and 2
  group_by(site_code, site_name, survey_date, date_time_survey, depth, method, species_name, common_name, size_class) %>%
  summarise(survey_total = sum(total)) %>% # sum blocks 1 and 2
  ungroup() %>%
  mutate(size_class = as.numeric(size_class),
         # correct for rectangle area
         survey_den = case_when(method == 1 ~ survey_total/500,
                                method == 2 ~ survey_total/100)) %>%
  as.data.frame() %>%
  select(site_code, site_name, survey_date, date_time_survey, depth, method, species_name, size_class, survey_total)

# merge the three!
full_df <- rbind(rls, kelp_rls, new_rls)
