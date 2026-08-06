#Figures for EDS report
# Load packages-----------------------------------------------------------------

library(readr)
library(tidyr)
library(dplyr) # for general data wrangling
library(lubridate) # for looking at dates
library(janitor) #to clean data frames
library(stringr) #to clean data frames
library(forcats) #to order factors

library(glmmTMB) # The swiss army knife of modeling packages
library(DHARMa) # inspect model residuals/check assumptions
library(ggeffects) # for extracting predictions and running post hoc tests
library(broom) # to run nested linear models
library(performance) #to check assumptions of my regressions on normal data
library(insight) #to "fix" performance package
library(lme4)

library(ggplot2) # Plotting data
library(patchwork) # Arrange multiple plots together
library(visreg) # plot model predictions
library(viridisLite) # colour palettes
library(metaDigitise) #to pull data from figures without raw values
library(rphylopic) #to put nice images on our plots
library(patchwork) #to stitch plots together
library(paletteer) #colour palette 

library(ggspatial)
library(sf) #geospatial package for points, lines, and polygons: to map the min distance from an outplant site
library(geosphere) #to map the min distance from an outplant site

#Figure fluff ---------------------------------------------------------------------------------------

#Fig 1a)
#using the rphylopic database we will pull an ID for a stock image of a chiton to graph
abalone_uuid <- get_uuid(name = "abalones")
abalone_uuid

#then we will make a small data frame to house the image containing x and y values relevant to our graph of interest so we can add it to our plot
silhouette_df_abalone <- data.frame(
  img_x = 0.07, 
  img_y = "helby_ne", 
  uuid = abalone_uuid)

#Fig 1b)
#using the rphylopic database we will pull an ID for a stock image of a chiton to graph
abalone_uuid_2 <- get_uuid(name = "abalones")
abalone_uuid_2

#then we will make a small data frame to house the image containing x and y values relevant to our graph of interest so we can add it to our plot
silhouette_df_abalone_2 <- data.frame(
  img_x = 0.06, 
  img_y = "execution_rock", 
  uuid = abalone_uuid_2)

#Fig 1a) --------------------------------------------------------------------------------------------

#the first figure is going to be plotting the change in density by status of site (outplant or not)
#to get change in density (i.e., through time) we can only really use the Jane Watson data set which has data points at 4 different time periods (before, before, during, after outplant). 
#we can then try to add the RLS data to this which has data from 2021-2025 but we must make assumptions on how compareable these are via methods, and we must lump close by sites or discard sites that aren't surveyed by both methods

#first make df have site level outplant information
density_outplant_v_non <- read.csv("data-processed/abalone_density_joined.csv") %>%
  mutate(outplant = case_when(str_detect(site, "scotts_bay") ~"yes", #add a column for outplant level
                              str_detect(site, "aguilar_point") ~"yes",
                              str_detect(site, "grappler") ~"yes",
                              str_detect(site, "ed_king_se") ~"yes",
                              str_detect(site, "sandford_sw") ~"yes",
                              str_detect(site, "helby_sw") ~"yes",
                              TRUE ~ "no")) #everything else = not outplant

#next, filter for just Jane data
Fig1_JW <- density_outplant_v_non%>%
  filter(method == "par_quad") %>% #clean this up later, but for now we know each study has a unique method so we can just filter for that 
  filter(site != "ed_king_?") #remove sites that dont record density in certain time periods

#lets plot density through time per site, we will get a linear relationship from this and extract the slope to get a change in density through time across a site
ggplot(data = Fig1_JW, aes(x = as.numeric(year), y = site_mean_density, colour = outplant)) +
  geom_point()+
  geom_smooth(method = "lm") +
  labs(x = "Year", y = "Abalone density (" ~m^-2*")", colour = "Recieved outplants?") + #add titles
  theme_classic()+ #white background
  guides(color = guide_legend(reverse = TRUE))+ #put yes status on the top of the legend
  scale_y_continuous(breaks = c(0,0.5,1,1.5,2,2.5, 3))+ #play with this so that we can see the full range of CI but also have real and relevant y axis limits (without we plot into -ve)
  facet_wrap(~site) + #show by site
  coord_cartesian(ylim = c(0, 3)) #use this to show relevant CI range

#now we need to pull the slope value to get an average change in abalone density group

#okay, lets now get the slope of each line in this plot to then plot change in desnity through time
slopes_Fig1_JW_df <- Fig1_JW %>%
  nest_by(site)%>% #like group by but gives it a key
  mutate(fit = list(lm(site_mean_density ~ year, data = data)))%>% #this should run multiple linear models on the seperate sites so we get a slope per site
  reframe(tidy(fit)) %>%
  filter(term == "year")%>%
  mutate(outplant = case_when(str_detect(site, "scotts_bay") ~"yes", #add a column for outplant level
                              str_detect(site, "aguilar_point") ~"yes",
                              str_detect(site, "grappler") ~"yes",
                              str_detect(site, "ed_king_se") ~"yes",
                              str_detect(site, "sandford_sw") ~"yes",
                              str_detect(site, "helby_sw") ~"yes",
                              TRUE ~ "no"))
  

slopes_Fig1_JW <- slopes_Fig1_JW_df %>%
  mutate(site = fct_reorder(site, estimate)) %>% #this orders the sites by largest estimate to smallest estimate
  ggplot(aes(x = estimate, y = site, colour = outplant)) +
  geom_point(size = 2)+
  geom_errorbar(
    aes(
      xmin = estimate - std.error,
      xmax = estimate + std.error,
      y = site), size = 1.2, width = 0.2)+
  theme_classic()+
  scale_color_manual(values = c("#140E3AFF", "#CD64B5FF"))+
  labs(x = "Rate of recovery", y = "Site", colour = "Recieved outplants?")+
  scale_y_discrete(labels = \(x) str_to_title(str_replace_all(x, "_", " ")))+
  geom_vline(xintercept = 0, linetype = "dotted") + #insert 0 line to show direction of slopes
  guides(color = guide_legend(reverse = TRUE)) +#put yes status on the top of the legend
  geom_phylopic(data= silhouette_df_abalone, aes(x=img_x, y = img_y, uuid = abalone_uuid), height = 3, inherit.aes = FALSE) #use the rphylopic package alongside the data frame we created to place the image in space (altering image size with the height function)

#plot
slopes_Fig1_JW

#okay when brain dead clean these names and axis titles

###model--------------------------------------------------------------------------------------------

#compare the two outplant groups mean density (not accounting for years?) -> this doesn't really explain a nice story... I think we need something that accounts for years
t.test(site_mean_density ~ outplant,
       data = Fig1_JW,
       alternative = "two.sided")
#p = 0.1041 no evidence of 

#I think I need to be doing lms 
Fig1_model <- lm(estimate ~ outplant, data= slopes_Fig1_JW_df)
#before we had site and year in here but year is inherently in the estimate value so I dont think that should be included
check_model(Fig1_model) #non normal 
summary(Fig1_model)

#or maybe 
Fig1_model <-lm(site_mean_density ~ outplant * year, data= Fig1_JW)
#interaction as outplant sites should vary dependant on the year in this data set (year since outplant being 2022) but this leaves a very non-normal check model fit... Additive does not because it shows that year and outplant are very co linear?
check_model(Fig1_model)
summary(Fig1_model)
#why significantly -ve relationship? is it because its intercept its predicting is so high based on the heavily weighted zero data???


#figure out how to visualize the check_model. Since the correlation is -1 the random effect structure is too complex for data
Fig1_model <-lmer(site_mean_density ~ outplant + year + ( year | site), data= Fig1_JW)
#this is the model we think
#interaction as outplant sites should vary dependant on the year in this data set (year since outplant being 2022) but this leaves a very non-normal check model fit... Additive does not because it shows that year and outplant are very co linear?
check_model(Fig1_model, show_dots = FALSE) # will not let me check for normalicy 
summary(Fig1_model)
#Fig 1b) -------------------------------------------------------------------------------------------

#Okay lets try this with the addition of RLS (Not sure how i feel about combining methods)
#first we need to find = sites and drop sites that cannot be considered close enough to be =
unique(Fig1_JW$site)
#list of sites that we have similar RLS sites to = ed_king_nw, execution_rock, helby_ne (questionable with Wizard N), scotts_bay, aguilar_point, ed_king_se, village_bay, wizard, grappler
#so we need to make these names the same, easiest way to do that is change the JW names because I already have code for changing those names in the EDS_all_data script

Fig1_RLS <- density_outplant_v_non%>% #make df just RLS so we can easily find unique and relevant site names to match with JW dataset
  filter(method == "par_belt") 

unique(Fig1_RLS$site) #view names to change

Fig1_JW_RLS <- density_outplant_v_non %>%
  filter(method == "par_quad" | method == "par_belt") %>% #this is filtering for JW and RLS by the methods they use
  mutate(site = case_when(str_detect(site, "Kii xin") ~"execution_rock", #rename close sites to match
                              str_detect(site, "Ed King SW Pyramid") ~"ed_king_nw",
                              str_detect(site, "Wizard Island North") ~"helby_ne", #questionable
                              str_detect(site, "Dodger Channel") ~"ed_king_se", #could make argument for taylor rock instead but this one faces the same swell direction
                              str_detect(site, "Wizard Island South") ~"wizard",
                              str_detect(site, "Kirby") ~"village_bay", #note the kirby in JW data set is actually on the other side of the island and this village site is closer to the RLS kirby
                              TRUE ~ site))%>%
  group_by(site)%>%
  filter(all(c(1988, 2023, 2024) %in% year)) %>% #filter out non JW and RLS sites based on year covered
  ungroup()%>%
  mutate(outplant = case_when(str_detect(site, "ed_king_se") ~"yes",
         TRUE ~ outplant)) #since we lumped the previous RLS site to this site we must now say that it is an outplant site by proxy of how close it is 

#Now lets plot density through time per site for the combined data set
ggplot(data = Fig1_JW_RLS, aes(x = as.numeric(year), y = site_mean_density, colour = outplant)) +
  geom_point()+
  geom_smooth(method = "lm") +
  labs(x = "Year", y = "Abalone density (" ~m^-2*")", colour = "Recieved outplants?") + #add titles
  theme_classic()+ #white background
  guides(color = guide_legend(reverse = TRUE))+ #put yes status on the top of the legend
  scale_y_continuous(limits = c(-1, 3), breaks = c(0,0.5,1,1.5,2,2.5, 3))+ #play with this so that we can see the full range of CI but also have real and relevant y axis limits (without we plot into -ve)
  facet_wrap(~site) #show by site

#Lets get the slopes of the lines again
slopes_Fig1_JW_RLS_df <- Fig1_JW_RLS %>%
  nest_by(site)%>% #like group by but gives it a key
  mutate(fit = list(lm(site_mean_density ~ year, data = data)))%>% #this should run multiple linear models on the seperate sites so we get a slope per site
  reframe(tidy(fit)) %>%
  filter(term == "year")%>%
  mutate(outplant = case_when(str_detect(site, "scotts_bay") ~"yes", #add a column for outplant level
                              str_detect(site, "aguilar_point") ~"yes",
                              str_detect(site, "grappler") ~"yes",
                              str_detect(site, "ed_king_se") ~"yes",
                              str_detect(site, "sandford_sw") ~"yes",
                              str_detect(site, "helby_sw") ~"yes",
                              TRUE ~ "no"))


slopes_Fig1_JW_RLS <- slopes_Fig1_JW_RLS_df %>%
  mutate(site = fct_reorder(site, estimate)) %>% #this orders the sites by largest estimate to smallest estimate
  ggplot(aes(x = estimate, y = site, colour = outplant)) +
  geom_point(size = 2.5)+
  geom_errorbar(
    aes(
      xmin = estimate - std.error,
      xmax = estimate + std.error,
      y = site), size = 1.2, width = 0.2)+
  theme_classic()+
  scale_color_manual(values = c("#140E3AFF", "#CD64B5FF"))+
  labs(x = "Rate of recovery", y = "Site", colour = "Recieved outplants?")+
  theme(axis.title = element_text(size = 18))+
  scale_y_discrete(labels = \(x) str_to_title(str_replace_all(x, "_", " ")))+
  geom_vline(xintercept = 0, linetype = "dotted") + #insert 0 line to show direction of slopes
  guides(color = guide_legend(reverse = TRUE)) +#put yes status on the top of the legend
  geom_phylopic(data= silhouette_df_abalone_2, aes(x=img_x, y = img_y, uuid = abalone_uuid_2), height = 1.5, inherit.aes = FALSE) #use the rphylopic package alongside the data frame we created to place the image in space (altering image size with the height function)

#plot
slopes_Fig1_JW_RLS

#fix names and axis titles
#use patchwork to stitch
#okay, here is where we use patchwork to combine plots
slopes_Fig1_JW / slopes_Fig1_JW_RLS +
  plot_annotation(tag_levels = list(c('(a)', '(b)')))
#maybe in the future I can make 1 legend and picture that captures this whole plot 

#or maybe more explicitly?
slopes_Fig1_JW / slopes_Fig1_JW_RLS +
  plot_annotation(tag_levels = list(c('(conservative)', '(combined)')))
#Fig 2a) -----------------------------------------------------------------------------

#here I want to plot all the relevant sites distance from one of 6 potential outplant sites

#read in relevant coordinate data
Fig2_coords <- read.csv("data-raw/Barkley_Sound_Sites_Decimal_Degrees_EDS_sites.csv") %>%
  clean_names()%>%
  mutate( #tidy site names
    site = case_when(
      row_number() == 1 ~ "ship_islands",
      row_number() == 2 ~ "ed_king_sw", #call relevant to outplant location
      row_number() == 3 ~ "ed_king_nw", #call relevant to outplant location
      row_number() == 4 ~ "seppings",
      row_number() == 5 ~ "cape_beale",
      row_number() == 6 ~ "lawton_point",
      row_number() == 7 ~ "whittlestone",
      row_number() == 8 ~ "execution_rock",
      row_number() == 9 ~ "self_point", #this sites gps location is near effingham inlet?
      row_number() == 10 ~ "helby_ne",  #call relevant to outplant location
      row_number() == 11 ~ "scotts_bay",
      row_number() == 12 ~ "aguilar_point",
      row_number() == 13 ~ "cia_rock",
      row_number() == 14 ~ "ed_king_se", #for the purposes of this project I am renaming this site to align with the outplant location name (it is the closest site in this study)
      row_number() == 15 ~ "kirby_point",
      row_number() == 16 ~ "village_bay",
      row_number() == 17 ~ "blowhole",
      row_number() == 18 ~ "prasiola",
      row_number() == 19 ~ "wizard",
      row_number() == 20 ~ "grappler",
      row_number() == 21 ~ "ed_king_?", #seems like the same place as ed_king_nw
      TRUE ~ site
    )) %>%
  filter(site != "ed_king_?") %>% #remove random site
  filter(site != "self_point") #remove self point as coord is incorrect and too late to contact Jane 
#NOTE I cannot locate a specific sanford island outplant location so for the purpose of this report I am going to neglect that data. Im not sure it will actually matter as all sites might be closer to another outplant location and


# Use package sf to find distance to closest outplant site - Convert to an sf object
sites_sf <- st_as_sf(
  Fig2_coords,
  coords = c("longitude_decimal_degrees", "latitude_decimal_degrees"),
  crs = 4326
)

# Split into outplant and non-outplant sites
outplant_sites <- sites_sf %>%
  filter(outplant == "yes")

non_outplant_sites <- sites_sf %>%
  filter(outplant == "no")

dist_mat <- st_distance(non_outplant_sites, outplant_sites) #make a distance matrix

non_outplant_sites$nearest_distance_m <-
  apply(dist_mat, 1, min) #find minimum distance from outplant

non_outplant_sites$nearest_outplant <-
  outplant_sites$site[apply(dist_mat, 1, which.min)]  #add which specific site it is

non_outplant_sites <- non_outplant_sites %>% #lets change this to km so its easier to graph and interpret
  mutate(
    nearest_distance_km = as.numeric(nearest_distance_m) / 1000
  )

#now I need to add nearest distance in km to the slope estimate data frame 
min_distance_join <- left_join(slopes_Fig1_JW_df, non_outplant_sites, by = c("site", "outplant"))%>%
  #filter(outplant != "yes") %>% #remove all sites that are outplants as we dont have a distance from relationship here
  filter(site != "self_point") %>% #have to remove again because joined
#Review self_point when have more time
  mutate(nearest_distance_km = case_when(is.na(nearest_distance_km) ~ 0, #add outplant sites
                                        TRUE ~ nearest_distance_km)) %>%
  mutate(nearest_outplant = case_when(str_detect(site, "scotts_bay") ~"scotts_bay", #include outplant locations as nearest outplants instead of NA's
                                      str_detect(site, "aguilar_point") ~"aguilar_point",
                                      str_detect(site, "grappler") ~"grappler",
                                      str_detect(site, "helby_sw") ~"helby_sw", #not in data set
                                      str_detect(site, "ed_king_se") ~"ed_king_se",
                                      TRUE ~ nearest_outplant))%>%
  filter(site != "aguilar_point")%>% #only one point and thats aguilar so no trend to detect
  filter(site != "grappler") #only one point and thats grappler so no trend to detect

#lets try to plot this now!
distance_from_outplant_Fig2_JW <- min_distance_join %>%
  ggplot(aes(x = nearest_distance_km, y = estimate, colour = nearest_outplant)) +
  geom_point() +
  #geom_smooth(aes(fill = nearest_outplant), method = "lm") +
  geom_smooth(method = "lm")+
  geom_errorbar(
    aes(
      ymin = estimate - std.error,
      ymax = estimate + std.error,
      x = nearest_distance_km),width = 0.2)+
  theme_classic()+
  labs(x = "Distance from closest large outplant (km)", y = "Rate of density change", colour = "Outplant location")+
  geom_hline(yintercept = 0, linetype = "dotted") +  #insert 0 line to show direction of slopes
  coord_cartesian(ylim = c(-0.07, 0.09)) #use this to show relevant error between points
#no scotts bay error becuase there are just two points
  

distance_from_outplant_Fig2_JW
#perhaps there are more applicable further away outplant sites for the sites labelled by ed_king_se. Important to remember that only around 10,000 abs outplanted here compared to the millions at scotts_bay, helby, and grappler

###model 2a)----------------------------------------------------------------------------------------
Fig2_model <- lm(estimate ~ nearest_distance_km + nearest_outplant, data= min_distance_join)
summary(Fig2_model)
#hmm this is predicting our slopes to all be positive which isn't the reality of this data set
#it probably does this because the weight of the data towards a site that was "outplanted" at a very different level is seen as an equal factor in the model to scotts and helby

#Fig 2b)--------------------------------------------------------------------------------------------
#HERE WE REMOVE THE ED KING OUTPLANT BECAUSE IT OUTPLANTED ADULTS ON THE MAGNITUDE OF 10,000 WHERE AS ITS OTHER TWO COMPARISONS OUTPLANTED MAINLY LARVAE AND JUVENILES ON THE SCALES OF MILLIONS. THEY ARE NOT THE SAME
b_Fig2_coords <- Fig2_coords%>%
  filter(site != "ed_king_se")

# Use package sf to find distance to closest outplant site - Convert to an sf object
b_sites_sf <- st_as_sf(
  b_Fig2_coords,
  coords = c("longitude_decimal_degrees", "latitude_decimal_degrees"),
  crs = 4326
)

# Split into outplant and non-outplant sites
b_outplant_sites <- b_sites_sf %>%
  filter(outplant == "yes")

b_non_outplant_sites <- b_sites_sf %>%
  filter(outplant == "no")

b_dist_mat <- st_distance(b_non_outplant_sites, b_outplant_sites) #make a distance matrix

b_non_outplant_sites$nearest_distance_m <-
  apply(b_dist_mat, 1, min) #find minimum distance from outplant

b_non_outplant_sites$nearest_outplant <-
  b_outplant_sites$site[apply(b_dist_mat, 1, which.min)]  #add which specific site it is

b_non_outplant_sites <- b_non_outplant_sites %>% #lets change this to km so its easier to graph and interpret
  mutate(
    nearest_distance_km = as.numeric(nearest_distance_m) / 1000
  )

#now I need to add nearest distance in km to the slope estimate data frame 
b_min_distance_join <- left_join(slopes_Fig1_JW_df, b_non_outplant_sites, by = c("site", "outplant"))%>%
  #filter(outplant != "yes") %>% #remove all sites that are outplants as we dont have a distance from relationship here
  filter(site != "self_point") %>% #have to remove again because joined
  #Review self_point when have more time
  mutate(nearest_distance_km = case_when(is.na(nearest_distance_km) ~ 0, #add outplant sites
                                         TRUE ~ nearest_distance_km)) %>%
  mutate(nearest_outplant = case_when(str_detect(site, "scotts_bay") ~"scotts_bay", #include outplant locations as nearest outplants instead of NA's
                                      str_detect(site, "aguilar_point") ~"aguilar_point",
                                      str_detect(site, "grappler") ~"grappler",
                                      str_detect(site, "helby_sw") ~"helby_sw", #not in data set
                                      str_detect(site, "ed_king_se") ~"ed_king_se",
                                      TRUE ~ nearest_outplant))%>%
  filter(site != "aguilar_point")%>% #only one point and thats aguilar so no trend to detect
  filter(site != "grappler") %>% #only one point and thats grappler so no trend to detect
  filter(site != "ed_king_se")

#lets try to plot this now!
b_distance_from_outplant_Fig2_JW <- b_min_distance_join %>%
  ggplot(aes(x = nearest_distance_km, y = estimate, colour = nearest_outplant)) +
  geom_point() +
  #geom_smooth(aes(fill = nearest_outplant), method = "lm") +
  geom_smooth(method = "lm")+
  geom_errorbar(
    aes(
      ymin = estimate - std.error,
      ymax = estimate + std.error,
      x = nearest_distance_km),width = 0.2)+
  theme_classic()+
  guides(color = guide_legend(reverse = TRUE))+ #put yes status on the top of the legend
  labs(x = "Distance from closest outplant (km)", y = "Rate of density change", colour = "Outplant location")+
  geom_hline(yintercept = 0, linetype = "dotted") +  #insert 0 line to show direction of slopes
  coord_cartesian(ylim = c(-0.07, 0.09)) #use this to show relevant CI range


b_distance_from_outplant_Fig2_JW

###model 2b)----------------------------------------------------------------------------------------
#lets try this again 
b_Fig2_model <- lm(estimate ~ nearest_distance_km + nearest_outplant, data= b_min_distance_join)
check_model(b_Fig2_model)
summary(b_Fig2_model)
#for every 1 km from the outplant site the rate of change decreases by -0.005 ab per year.
#there is moderate evidence to suggest that for every 1 km from the outplant site, scotts bay has a more positve rate of change.
#Based off this evidence, particularly the scotts bay site, say unit at witch crosses -ve slope. 
#explain this in discussion based off ecological relevance 

#In the future--------------------------------------------------------------------------------------
#Hmm, could we instead do this with mean site abalone density using data sets post outplanting (2021-2026)
