#MSc figures 
#Oct 20th, 2026

#load packages
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
library(TMB)
library(visreg)
library(fitdistrplus)
library(emmeans)

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


#Figure fluff ----
#Fig 3
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

#Fig 2
#using the rphylopic database we will pull an ID for a stock image of a chiton to graph
abalone_uuid_3 <- get_uuid(name = "abalones")
abalone_uuid_3

#then we will make a small data frame to house the image containing x and y values relevant to our graph of interest so we can add it to our plot
silhouette_df_abalone_3 <- data.frame(
  img_x = 6, 
  img_y = 0.06, 
  uuid = abalone_uuid_3)


#Fig 1: Map ----
## JW + annual ----
JW_coords <- read.csv("data-raw/Barkley_Sound_Sites_Decimal_Degrees_EDS_sites.csv")
#must do this in order because of slice usage
JW_coords <- janitor::clean_names(JW_coords)%>% 
  mutate(site = str_to_lower(str_replace_all(site, " ", "_")), #tidy the site names
         site = str_replace_all(site, "'", ""))%>%
  filter(site != "self_point") %>% #remove random site until get gps
  slice(-20) #remove ed_king_? which only had data from one year

## KRI & Parks ----
KRI_coords <- read.csv("data-raw/KRI_coords.csv") 

KRI_coords <- janitor::clean_names(KRI_coords) %>%
  mutate(
    site = str_to_lower(site),
    site = str_replace_all(site, "\\s+", "_"), #changes all white spaces in the names 
    site = str_replace_all(site, "'", ""),
    outplant = "no")

## Storm Coast ----
SC_coords <- read.csv("data-raw/StormCoast_coords.csv") 

SC_coords <- janitor::clean_names(SC_coords) %>%
  rename(site = site_id)%>%
  mutate(
    site = str_to_lower(site),
    site = str_replace_all(site, "\\s+", "_"), #changes all white spaces in the names 
    site = str_replace_all(site, "'", ""),
    outplant = "no")

SC_coords_clean <- SC_coords %>%
  dplyr::select(site, transect_id, latitude_decimal_degrees, longitude_decimal_degrees)%>%
  mutate(outplant = case_when(str_detect(site, "eagle_bay") ~"yes", #add a column for outplant level
                              str_detect(site, "aguilar") ~"yes",
                              str_detect(transect_id, "GOT01") ~"yes",
                              str_detect(transect_id, "EKE02") ~"yes",
                              str_detect(transect_id, "HSW01") ~"yes",
                              TRUE ~ "no")) #everything else = not outplant
## RLS + KCCA ----  
#load RLS 2026 data 
RLS_coords_2026 <- read_csv("data-raw/RLS_2026_only_bonus.csv")

#clean
RLS_coords_2026 <- janitor::clean_names(RLS_coords_2026) %>%
  slice(-1) #remove first RLS row that is remnant RLS formatting

#identify unique site locations (effingham west and Baeria north island south side)
RLS_coords_clean <- RLS_coords_2026 %>%
  distinct(site_name, site_no, latitude, longitude) %>% #notice  even if two transects at a single site the lat, long stays same so dont need depth!
  rename(site = site_name)%>%
  mutate(
    site = str_to_lower(site),
    site = str_replace_all(site, "\\s+", "_"), #changes all white spaces in the names 
    site = str_replace_all(site, "'", ""),
    outplant = "no")

RLS_coords_2026 <- RLS_coords_clean %>%
  mutate(outplant = case_when(str_detect(site, "eagle_bay") ~"yes", #add a column for outplant level
                              str_detect(site, "aguilar_point") ~"yes",
                              str_detect(site, "goby_town") ~"yes",
                              str_detect(site, "EKE02") ~"yes",
                              str_detect(site, "black_fish_rock") ~"yes",
                              TRUE ~ "no")) #everything else = not outplant

#load RLS data 2021-2025 which contains two sites not included in the above
RLS_coords_2025 <- read_csv("data-raw/RLS_2025.csv")

#clean
RLS_coords_2025 <- janitor::clean_names(RLS_coords_2025) 

#identify unique site locations (effingham west and Baeria north island south side)
RLS_coords_clean_2025 <- RLS_coords_2025 %>%
  distinct(site_name, site_code, latitude, longitude) %>% #notice  even if two transects at a single site the lat, long stays same so dont need depth!
  rename(site = site_name)%>%
  rename(site_no = site_code)%>%
  mutate(
    site = str_to_lower(site),
    site = str_replace_all(site, "\\s+", "_"), #changes all white spaces in the names 
    site = str_replace_all(site, "'", ""),
    outplant = "no")%>%
  filter(site == "effingham_west" | site == "baeria_rocks_north_island_southside") 
         
RLS_coords <- bind_rows(RLS_coords_clean_2025, RLS_coords_2026)%>%
  dplyr::select(site, latitude, longitude, outplant)

#add KCCA RlS files
KCCA_coords <- read_csv("data-raw/KCCA_RLS_all.csv")

#clean
KCCA_coords <- janitor::clean_names(KCCA_coords) 

#identify unique site locations (effingham west and Baeria north island south side)
KCCA_coords <- KCCA_coords %>%
  distinct(site_name, latitude, longitude) %>%
  rename(site = site_name)%>%
  mutate(
    site = str_to_lower(site),
    site = str_replace_all(site, "\\s+", "_"), #changes all white spaces in the names 
    site = str_replace_all(site, "'", ""),
    outplant = "no") %>%
  filter(site != "sand_town")


## DFO ----
#load data
DFO_coords <- read.csv("data-raw/DFOKelp2021-2025.csv")
# clean
DFO_coords <- janitor::clean_names(DFO_coords)

DFO_coords <- DFO_coords%>% 
  rename(site = sites) %>%
  mutate(site = str_to_lower(str_replace_all(site, " ", "_")), #tidy the site names
         site = str_replace_all(site, "'", "")) 
# CHECK: NEED TO ASK SANDIE FOR GPS OF NEW SITES AND GENERAL RUN DOWN OF SURVEY AREA ETC
#Fig 3: JW rate ----

#first make df have site level outplant information
rate <- read.csv("data-processed/all_abalone_density.csv") %>%
  mutate(outplant = case_when(str_detect(site, "scotts_bay") ~"yes", #add a column for outplant level
                              str_detect(site, "aguilar_point") ~"yes",
                              str_detect(site, "grappler") ~"yes",
                              str_detect(site, "ed_king_se") ~"yes",
                              str_detect(site, "sandford_sw") ~"yes",
                              str_detect(site, "helby_sw") ~"yes",
                              TRUE ~ "no")) %>% #everything else = not outplant
  mutate(site = case_when(str_detect(site, "_se") ~ str_replace(site, "_se", "_southeast"), #rename sites that are just directions to their full name
                          str_detect(site, "_ne") ~ str_replace(site, "_ne", "_northeast"),
                          str_detect(site, "_nw") ~ str_replace(site, "_nw", "_northwest"),
                          str_detect(site, "_sw") ~ str_replace(site, "_sw", "_southwest"),
                          TRUE ~ site))

#filter for just JW data
Fig3_rate <- rate%>%
  filter(source == "Jane_Watson") %>% #filter for rate data (1988 - 2026)
  filter(site != "ed_king_?") %>% #remove sites that dont record density in certain time periods
  mutate(site = factor(site,           #manually (ew) change names to read largest positive slope to largest negative slopes 
                       levels = c(
                         "scotts_bay",
                         "grappler",
                         "taylor",
                         "aguilar_point",
                         "execution_rock",
                         "wizard",
                         "helby_northeast",
                         "self_point",
                         "whittlestone",
                         "prasiola",
                         "cape_beale",
                         "ed_king_northwest",
                         "ship_islands",
                         "seppings",
                         "lawton_point",
                         "cia_rock",
                         "kirby_point",
                         "village_bay",
                         "blowhole",
                         "ed_king_southeast",
                         "ed_king_southwest")))


#lets plot density through time per site, we will get a linear relationship from this and extract the slope to get a change in density through time across a site
Fig3_random_sites <- Fig3_rate %>%
  filter(!method == "perp_belt") %>%
  ggplot(aes(x = as.numeric(year), y = site_mean_density, colour = outplant)) +
  geom_point()+
  geom_smooth(method = "lm") +
  labs(x = "Year", y = "Abalone density (" ~m^-2*")", colour = "Recieved outplants?") + #add titles
  theme_classic()+ #white background
  theme(axis.title = element_text(size = 20), #increase axis text size
        legend.title = element_text(size = 18), #increase legend text size
        legend.text = element_text(size = 16))+
  scale_color_manual(values = c("grey65", "#CD64B5FF"))+
  guides(color = guide_legend(reverse = TRUE))+ #put yes status on the top of the legend
  scale_y_continuous(breaks = c(0,0.5,1,1.5,2,2.5, 3))+ #play with this so that we can see the full range of CI but also have real and relevant y axis limits (without we plot into -ve)+
  scale_x_continuous(breaks = c(1988,2007,2022)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ #adjust angle to fit without overlap
  facet_wrap(~site, #show by site
             ncol = 5, #to get 7 plots across the 3 rows
             labeller = labeller(site = \(x) str_to_title(str_replace_all(x, "_", " ")))) + # edit names) 
  theme(strip.text = element_text(size = 6.5)) +# Decrease size value
  coord_cartesian(ylim = c(0, 3)) #use this to show relevant CI range

#plot 
Fig3_random_sites

#pull in the annual sites
Fig3_annual_sites <- Fig3_rate %>%
  filter(!method == "par_quad") %>%
  ggplot(aes(x = as.numeric(year), y = site_mean_density, colour = outplant)) +
  geom_point(shape = 0)+
  geom_smooth(method = "lm") +
  labs(x = "Year", colour = "Recieved outplants?") + #add titles
  theme_classic()+ #white background
  theme(axis.title = element_text(size = 20), #increase axis text size
        axis.title.y = element_blank()) + #remove y title for combined plot
  scale_color_manual(values = c("grey65", "#CD64B5FF"))+
  scale_y_continuous(breaks = c(0,0.5,1,1.5,2,2.5, 3))+ #play with this so that we can see the full range of CI but also have real and relevant y axis limits (without we plot into -ve)+
  scale_x_continuous(breaks = c(1987,2007,2026)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ #adjust angle to fit without overlap
  facet_wrap(~site, #show by site
             ncol = 2, #to get 7 plots across the 3 rows
             labeller = labeller(site = \(x) str_to_title(str_replace_all(x, "_", " ")))) + # edit names) 
  theme(strip.text = element_text(size = 10)) + # Decrease size value
  coord_cartesian(ylim = c(0, 3)) #use this to show relevant CI range

#plot 
Fig3_annual_sites

#combine figs
Fig3_combined <- Fig3_random_sites | Fig3_annual_sites +
  plot_layout(widths = c(3, 1.5)) &
  theme(legend.position = "none")

Fig3_combined
#ggsave("figures/MSc/Fig3a.png", device = "png", height = 9, width = 12, dpi = 400)
#commented out so that when you run the code it doesn't save whats in your plot accidentally :)

#okay, lets now get the slope of each line in this plot to then plot change in density through time
slopes_Fig3 <- Fig3_rate %>%
  mutate(site = case_when(site == "wizard" & method == "perp_belt" ~ "wizard_north",
    TRUE ~ site)) %>%
  nest_by(site)%>% #like group by but gives it a key
  mutate(fit = list(lm(site_mean_density ~ year, data = data)))%>% #this should run multiple linear models on the separate sites so we get a slope per site
  reframe(tidy(fit)) %>%
  filter(term == "year")%>%
  mutate(outplant = case_when(str_detect(site, "scotts_bay") ~"yes", #add a column for outplant level
                              str_detect(site, "aguilar_point") ~"yes",
                              str_detect(site, "grappler") ~"yes",
                              str_detect(site, "ed_king_se") ~"yes",
                              str_detect(site, "sandford_sw") ~"yes",
                              str_detect(site, "helby_sw") ~"yes",
                              TRUE ~ "no")) %>%
  mutate(method = case_when(str_detect(site, "wizard_north") ~"perp_belt", #add a column for method level
                              str_detect(site, "taylor") ~"perp_belt",
                              TRUE ~ "par_quad"))
  

slopes_Fig3_plot <- slopes_Fig3 %>%
  mutate(site = fct_reorder(site, estimate)) %>% #this orders the sites by largest estimate to smallest estimate
  ggplot(aes(x = estimate, y = site, colour = outplant, shape = method)) +
  geom_point(size = 3)+
  scale_shape_manual(values = c(16, 0)) +
  geom_errorbar(
    aes(
      xmin = estimate - std.error,
      xmax = estimate + std.error,
      y = site), size = 1.2, width = 0.2)+
  theme_classic()+
  theme(axis.title = element_text(size = 19), #increase axis text size
        legend.title = element_text(size = 17), #increase legend text size
        legend.text = element_text(size = 15),
        legend.position = "none") + #remove legend up so can have one legend in patched figure 1
  scale_color_manual(values = c("grey65", "#CD64B5FF"))+
  labs(x = "Rate of abalone density change through time", y = "Site", colour = "Recieved outplants?")+
  scale_y_discrete(labels = \(x) str_to_title(str_replace_all(x, "_", " ")))+
  scale_x_continuous(breaks = c(-0.025, 0.00, 0.05)) + #add -ve x landmark 
  geom_vline(xintercept = 0, linetype = "dotted") + #insert 0 line to show direction of slopes
  #scale_linetype_manual(name = "Recovery?",values = c("Recovering?" = "dotted"))+ # try to add legend item for dashed line - not working
  guides(color = guide_legend(reverse = TRUE)) #put yes status on the top of the legend
#geom_phylopic(data= silhouette_df_abalone, aes(x=img_x, y = img_y, uuid = abalone_uuid), height = 4, inherit.aes = FALSE) #use the rphylopic package alongside the data frame we created to place the image in space (altering image size with the height function)

#plot
slopes_Fig3_plot

#bring all Fig 3 together
Fig3_combined / slopes_Fig3_plot +
  plot_layout(guides = "collect")+ #bring legend to middle of the patchwork grid
  plot_annotation(tag_levels = list(c('(a)', '(b)', '(c)')))

#model
Fig3_model_test <-lmer(site_mean_density ~ outplant + factor(year) + ( 1| site), data= Fig3_rate)
check_model(Fig3_model_test, show_dots = FALSE) # will not let me check for normalcy 
plot(simulateResiduals(Fig3_model_test))


Fig3_model <-lm(site_mean_density ~ outplant * year, data= Fig3_rate)
check_model(Fig3_model_test, show_dots = FALSE) # will not let me check for normalcy 
plot(simulateResiduals(Fig3_model_test))
summary(Fig3_model)


simulate_residuals(Fig3_model)

#model
descdist(Fig3_glm$site_mean_density) #looks like gamma but has 0s so try tweedie
Fig3_glm <- Fig3_rate %>%
  #filter(method != "perp_belt") %>% #remove annual data for modelling 
  #notice that if you add the filter above you take a lot of the weight out of the non_outplant data that is rich in the annual surveys and this drives a non-significant pattern with the interaction but generates one with year
  mutate(year_decade = (year - mean(year, na.rm = TRUE)) / 10) %>% #attempt to make year continuous and see "change per decade"
  mutate(pre_post = case_when(year > 2003 ~ "pre",
                          TRUE ~ "post"))
Fig3_model_GLM <- glmmTMB(site_mean_density ~ outplant * year_decade * pre_post + (1| site), family = tweedie(link = "log"), data = Fig3_glm) 
#tweedie handles the 0s in the dataset, and can handle the random factor and the unbalanced survey design
plot(simulateResiduals(Fig3_model_GLM))
summary(Fig3_model_GLM) #have to exp to interpret 
exp(0.66322) #back transform the estimate for outplant_yes = 1.941032 1.94x higher density at outplant sites
exp(0.25386) #back transform the estimate for decade = 1.288991 abs per m2 per decade growth at non-outplant sites

#test site as main effect for better slope estimates!
Fig3_model_GLM <- glmmTMB(site_mean_density ~  year_decade + site, family = tweedie(link = "log"), data = Fig3_glm)
plot(simulateResiduals(Fig3_model_GLM))
summary(Fig3_model_GLM) 
# Biggest question is don't I then want to calculate the slopes using this distribution not lm as done for the Fig 3C, so how do I do this?

#Isa model
Fig3_model_GLM <- glmmTMB(site_mean_density ~ outplant + year_decade + pre_post + outplant*year_decade + outplant*pre_post + year_decade*pre_post + (1| site), family = tweedie(link = "log"), data = Fig3_glm)
plot(simulateResiduals(Fig3_model_GLM))
summary(Fig3_model_GLM) #have to exp to interpret
