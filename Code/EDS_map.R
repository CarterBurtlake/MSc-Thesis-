#Map of EDS report sites
#August 2026
#Carter Burtlake

# Load packages ----------------------------------------------------------------
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
library(paletteer) #colour palette 

library(ggspatial)
library(sf)

# Load data --------------------------------------------------------------------

JW_coords <- read.csv("data-raw/Barkley_Sound_Sites_Decimal_Degrees_EDS_sites.csv") %>%
  janitor::clean_names()%>%
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

# transform coords into actually mappable coords
coords <- JW_coords %>%
  st_as_sf(coords = c("longitude_decimal_degrees", "latitude_decimal_degrees")) %>%
  st_set_crs(4326) %>%
  left_join((JW_coords %>% select(site, latitude_decimal_degrees, longitude_decimal_degrees)), by = "site")%>%
  rename(latitude = "latitude_decimal_degrees") %>%
  rename(longitude = "longitude_decimal_degrees") %>%
  mutate(
    site_label = stringr::str_to_title(gsub("_", " ", site))
  )

# Mapping ----------------------------------------------------------------------

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
                             size = 3,
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
    annotation_scale(location = "br", width_hint = 0.4) +
    annotation_north_arrow(location = "br", which_north = "true", 
                           pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering)
}


# Make a map ----
map(lat_min = -125.245, #set map boundary
          lat_max = -125.09, 
          long_min = 48.76, 
          long_max = 48.87, 
          coord_data = coords, #add data
         # nh4_var = mean_density, #use later for density
          point_size = 6, #set size intentionally large to obscure sites
          kelp_var = outplant, #label outplant sites
          map_file = hakai_map) #pick map file

#ggsave("figures/map.png", device = "png", height = 9, width = 12, dpi = 400)

