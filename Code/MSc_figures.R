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
  slice(-20) %>% #remove ed_king_? which only had data from one year
  mutate(source = "Jane_Watson") %>%
  rename(latitude = latitude_decimal_degrees,  #change column names to match other dfs
         longitude = longitude_decimal_degrees)
  
JW_coords <- JW_coords %>% #fix to match other dfs
  mutate(outplant = case_when(str_detect(site, "scotts_bay") ~"yes", 
                              #str_detect(site, "aguillar_point") ~"yes",      #not primary outplant
                              str_detect(site, "grappler") ~"yes",
                              #str_detect(site, "taylor_ed_king") ~"yes", #not primary outplant
                              TRUE ~ "no")) #everything else = not outplant
## KRI & Parks ----
KRI_coords <- read.csv("data-raw/KRI_coords.csv") 

KRI_coords <- janitor::clean_names(KRI_coords) %>%
  mutate(
    site = str_to_lower(site),
    site = str_replace_all(site, "\\s+", "_"), #changes all white spaces in the names 
    site = str_replace_all(site, "'", ""),
    outplant = "no") %>%
  mutate(source = "KRI")

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
                              #str_detect(site, "aguilar") ~"yes",      #not primary outplant
                              str_detect(transect_id, "GOT01") ~"yes",
                              #str_detect(transect_id, "EKE02") ~"yes", #not primary outplant
                              str_detect(transect_id, "HSW01") ~"yes",
                              TRUE ~ "no")) %>% #everything else = not outplant
  mutate(source = "Storm_Coast") %>%
  rename(latitude = latitude_decimal_degrees,
         longitude = longitude_decimal_degrees)

SC_coords <- SC_coords_clean[-2]
#CHECK: one site missing???

## RLS ----  
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
                              #str_detect(site, "aguilar_point") ~"yes", #not primary outplant
                              str_detect(site, "goby_town") ~"yes",
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
  dplyr::select(site, latitude, longitude, outplant) %>%
  mutate(source = "RLS")

## KCCA ----
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
  filter(site != "sand_town") %>%
  mutate(source = "KCCA")


## DFO ----
#load data
DFO_coords <- read.csv("data-raw/DFOKelp2021-2025.csv")
# clean
DFO_coords <- janitor::clean_names(DFO_coords)

DFO_coords <- DFO_coords%>% 
  rename(site = sites) %>%
  mutate(site = str_to_lower(str_replace_all(site, " ", "_")), #tidy the site names
         site = str_replace_all(site, "'", "")) %>%
  mutate(outplant = case_when(str_detect(site, "scotts_bay") ~"yes", #add a column for outplant level
                              str_detect(site, "grappler_mouth") ~"yes",
                              TRUE ~ "no")) %>% #everything else = not outplant
  mutate(source = "DFO") %>%
  rename(latitude = latitude_decimal_degrees,  #change column names to match other dfs
         longitude = longitude_decimal_degrees) %>%
  filter(site != "new_taylor_island_2022", #remove for now while dont have GPS!
         site != "new_flemming_112_2022",
         site != "tzartus_116") #filter this site out as replaced in 2022

# CHECK: NEED TO ASK SANDIE FOR GPS OF NEW SITES AND GENERAL RUN DOWN OF SURVEY AREA ETC

## Combined coords ----
coords <- rbind(JW_coords, DFO_coords, RLS_coords, KRI_coords, SC_coords, KCCA_coords)

### Map ----
# transform coords into actually mappable coords
coords_clean <- coords %>%
  st_as_sf(coords = c("longitude", "latitude"),
           remove = FALSE) %>%
  st_set_crs(4326)%>%
  mutate(
    site_label = stringr::str_to_title(gsub("_", " ", site)))

# Convert to projected coordinates (metres)
coords_projected <- coords_clean %>%
  st_transform(32610)

# Identify observations within 200 m
nearby_sites <- st_is_within_distance(
  coords_projected,
  coords_projected,
  dist = 200
)

nearby_pairs <- data.frame(
  pair = rep(seq_along(nearby_sites),
             lengths(nearby_sites)),
  match = unlist(nearby_sites)
) %>%
  filter(pair < match) %>%
  mutate(
    site_1 = coords$site[pair],
    site_2 = coords$site[match],
    source_1 = coords$source[pair],
    source_2 = coords$source[match]
  )
#okay lots are close... lets plot them all and then make adjustments from there
# Load GREAT shapefile 
hakai_map <- sf::st_read("data-raw/Hakai_coast/COAST_TEST2.shp") %>%
  st_sf() %>%
  st_set_crs(4326)

# honestly not 100% sure what this does
sf_use_s2(FALSE)

# Define colours
blue <- paste("#b9d1df", sep="")

# Making RLS maps -----
# make a map function
map <- function(lat_min, lat_max, long_min, long_max, 
                coord_data, nh4_var, kelp_var, point_size, map_file, 
                white_background = TRUE, invert = FALSE) {
  
  # invert land and sea colours if I use the potato map
  if(invert == FALSE){
    sea <- blue
    land <- "white"
  }
  if(invert == TRUE){
    sea <- "white"
    land <- blue
  }
  
  # specify black or white background
  if(white_background == TRUE){
    background <- "white"
    features <- "black"
  }
  if(white_background == FALSE){
    background <- "black"
    features <- "white"
  }
  
  ggplot() +
    geom_sf(data = map_file, fill = land, colour = sea) +
    # add points
    geom_sf(data = coord_data,
            alpha = 0.9,
            size = point_size,
            stroke = 2,
            aes(
              fill = {{nh4_var}}, #use this fill later for desnity scores. Also nice to not have big dots over the map as we want to obscure exact site location
              shape = {{kelp_var}},
              colour = {{kelp_var}}))  +
    viridis::scale_fill_viridis(option="magma", direction = -1,
                                guide = guide_colorbar(frame.colour = features, ticks.colour = features)) + #this will be used later to have a continuous colour fill for density scores
    coord_sf(xlim = c(lat_min, lat_max), ylim = c(long_min, long_max), expand = FALSE) + 
    ggrepel::geom_text_repel(data = {{coord_data}},
                             aes(x = longitude, y = latitude, 
                                 label = site_label), nudge_x = 0.001, #get name beside symbol
                             direction = "y",
                             hjust = 0,
                             size = 1,
                             segment.color = NA) +
    labs(fill = "Abalone density", pch = "Outplant site?") + #leave fill for later, for now keep to see map behind obscured data point. pch is for the shape
    scale_shape_manual(
      name = "Outplant site?",
      values = c("yes" = 21, "no" = 21), #star and circle
      breaks = c("yes", "no")
    ) +
    scale_colour_manual(
      values = c("yes" = "#CD64B5FF", "no" = "grey65"),
      breaks = c("yes", "no"),
      guide = "none"
    ) +
    guides(pch = guide_legend(override.aes = 
                                list(colour = c("#CD64B5FF", "grey65")))) +
    # Themes
    theme_bw() +
    theme(
      # panel stuff
      panel.background = element_rect(fill = sea),
      panel.grid.major = element_line(color = sea),
      panel.border = element_rect(fill = NA, colour = features),
      # remove axis
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.ticks.length = unit(0, "pt"),
      plot.title = NULL,
      plot.margin=grid::unit(c(0,0,0,0), "mm"),
      # Specify legend options
      legend.background = element_rect(color = NA, fill = background),  
      legend.key = element_rect(color = background,  fill = background),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = 24, color = features),  
      legend.title = element_text(size = 24, face = "bold", hjust = 0, color = features),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL,
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = 30, color = features),  
      strip.text.y = element_text(size = 30, color = features,angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = background, fill = background),  
    ) +
    annotation_scale(location = "br", width_hint = 0.3) +
    annotation_north_arrow(location = "br", which_north = "true", 
                           pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
}


# Make a map ----
map(lat_min = -125.4, #set map boundary
    lat_max = -124.99, 
    long_min = 48.73, 
    long_max = 49.00, 
    coord_data = coords_clean, #add data
    # nh4_var = mean_density, #use later for density
    point_size = 1, #set size intentionally large to obscure sites
    kelp_var = outplant, #label outplant sites
    map_file = hakai_map) #pick map file

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
  scale_x_continuous() +
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



# Nov 1 check jane analysis

Fig3_rate
str(Fig3_rate)
count(Fig3_rate, year)

Fig3_glm <- Fig3_rate %>%
  #filter(method != "perp_belt") %>% #remove annual data for modelling 
  #notice that if you add the filter above you take a lot of the weight out of the non_outplant data that is rich in the annual surveys and this drives a non-significant pattern with the interaction but generates one with year
  mutate(year_decade = (year - mean(year, na.rm = TRUE)) / 10) %>% #attempt to make year continuous and see "change per decade"
  mutate(pre_post = case_when(year > 2003 ~ "post",
                              TRUE ~ "pre")) %>%
  mutate(outplant = case_when(site == "ed_king_southeast" ~ "no",
                              site == "aguilar_point" ~ "no",
                              TRUE ~ outplant))

str(Fig3_glm)
JW_abalone_density_7year

rate <- read.csv("data-processed/all_abalone_density.csv")
str(rate)
levels(as.factor(rate$method))
count(rate,method,site)
## par_quad
## perp_belt

str(Fig3_rate)
JaneAnnual=Fig3_glm[Fig3_glm$method == "perp_belt", ]
str(JaneAnnual)

count(JaneAnnual,pre_post, year)

levels(as.factor(JaneAnnual$pre_post))
JaneAnnual$pre_post=as.factor(JaneAnnual$pre_post)
# Make sure Pre/Post is ordered correctly
JaneAnnual$pre_post <- factor(
  JaneAnnual$pre_post,
  levels = c("pre", "post"))

#--------------------------------------------------
# Calculate regression equations for each Site × Period
#--------------------------------------------------

Equation_Data <- JaneAnnual %>%
  group_by(site, pre_post) %>%
  summarise(
    model = list(lm(site_mean_density ~ year)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    intercept = coef(model)[1],
    slope = coef(model)[2],
    r2 = summary(model)$r.squared,
    
    equation = paste0(
      pre_post,
      ": y = ",
      round(intercept, 2),
      ifelse(slope >= 0, " + ", " - "),
      abs(round(slope, 2)),
      "x",
      ", R² = ",
      round(r2, 2)
    )
  ) %>%
  ungroup()
Equation_Data

#--------------------------------------------------
# Give equations different vertical positions
#--------------------------------------------------

Equation_Data <- Equation_Data %>%
  group_by(site) %>%
  mutate(
    label_y = c(
      max(JaneAnnual$site_mean_density, na.rm = TRUE),
      max(JaneAnnual$site_mean_density, na.rm = TRUE) * 0.90
    )
  ) %>%
  ungroup()


#--------------------------------------------------
# Plot
#--------------------------------------------------

ggplot(
  JaneAnnual,
  aes(
    x = year,
    y = site_mean_density,
    colour = pre_post)) +
  geom_point(
    size = 2.5,
    alpha = 0.8) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 1) +
  geom_text(
    data = Equation_Data,
    aes(
      x = -Inf,
      y = label_y,
      label = equation,
      colour = pre_post),
    inherit.aes = FALSE,
    hjust = -0.05,
    size = 4) +
  facet_wrap(
    ~ site,
    nrow = 1) +
  labs(
    x = "X variable",
    y = "Y variable",
    colour = NULL) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    strip.text = element_text(
      size = 11,
      face = "bold"),
    legend.position = "bottom")


### jane annual surveys

JaneRandom=Fig3_glm[Fig3_glm$method == "par_quad", ]
count(JaneRandom,site,year_decade)
str(JaneRandom)
hist(JaneRandom$site_mean_density)

JaneRandom_Model1 <- glmmTMB(site_mean_density ~ pre_post*outplant+ (1| site), family = tweedie(link = "log"), data = JaneRandom)
#can you bake in "in that pre category are outplant yes different than outplant no"
#make a boxplot of these 4 variables to determine this

summary(JaneRandom_Model1)
emm=emmeans(JaneRandom_Model1,~pre_post|outplant,type="response")
pairs(emm)


#to maintain degrees of freedom but just compare the effect of quadrat vs belt 
emmean_simple_pre_post <- emmeans(JaneRandom_Model1, pairwise ~ c(outplant, pre_post), simple = "pre_post")
emmean_simple_pre_post
#make a boxplot of these 4 variables

JaneRandom$year=as.factor(JaneRandom$year)
aggregate(site_mean_density~pre_post+outplant+year, data=JaneRandom,FUN=mean)

aggregate(site_mean_density~pre_post+outplant+year, data=JaneRandom,FUN=mean)

JaneRandom %>% 
  group_by(pre_post,outplant,year)%>%
summarise(
  mean_density=mean(site_mean_density),
  sd_density=sd(site_mean_density),
  .groups="drop")
  
  
JaneRandom %>% 
  group_by(year)%>%
  summarise(
    mean_density=mean(site_mean_density),
    sd_density=sd(site_mean_density),
    .groups="drop")

#check via boxplot with medium outplant sites (n=4)
#make df without annual but with medium outplant sites
medium_outplant_JW <- Fig3_rate %>%
  #filter(method != "perp_belt") %>% #remove annual data for modelling 
  #notice that if you add the filter above you take a lot of the weight out of the non_outplant data that is rich in the annual surveys and this drives a non-significant pattern with the interaction but generates one with year
  mutate(year_decade = (year - mean(year, na.rm = TRUE)) / 10) %>% #attempt to make year continuous and see "change per decade"
  mutate(pre_post = case_when(year > 2003 ~ "post",
                              TRUE ~ "pre")) %>%
  filter(method == "par_quad")
ggplot(medium_outplant_JW,
       aes(x = factor(year), y = site_mean_density, fill = outplant)) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(values = c("yes" = "lightblue",
                               "no" = "lightpink"), 
                    breaks = c("yes", "no")) +
  labs(
    x = "Year",
    y = expression("Abalone density (" ~m^-2*")"),
    fill = "Outplant status"
  ) +
  theme_minimal()

#check via boxplot with primary outplant sites (n=2)
ggplot(JaneRandom, aes(x = year, y = site_mean_density, fill = outplant)) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(values = c("yes" = "lightblue",
                               "no" = "lightpink"),
                    breaks = c("yes", "no")) +
  labs(
    x = "Year",
    y = expression("Abalone density (" ~m^-2*")"),
    fill = "Outplant status"
  ) +
  theme_minimal()
#based off this, it almost makes sense to not even consider pre / post as a grouping, and instead just have year and outplant


#model based off main outplants
JaneRandom_Model2 <- glmmTMB(site_mean_density ~ year*outplant+ (1| site), family = tweedie(link = "log"), data = JaneRandom)
plot(simulateResiduals(JaneRandom_Model2)) #looks good
summary(JaneRandom_Model2)

#model based off medium outplants
JaneRandom_Model3 <- glmmTMB(site_mean_density ~ as.factor(year)*outplant+ (1| site), family = tweedie(link = "log"), data = medium_outplant_JW)
plot(simulateResiduals(JaneRandom_Model3)) #looks good
summary(JaneRandom_Model3)
## Nov 1 Jane check end 










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
  
#test
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

#test
Fig3_model_GLM <- glmmTMB(site_mean_density ~ 0+pre_post + year_decade + outplant + (1| site), family = tweedie(link = "log"), data = Fig3_glm)

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
