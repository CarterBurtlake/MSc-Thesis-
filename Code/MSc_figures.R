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
library(igraph) #visualization


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
  filter(site != "self_point",  #remove random site until get gps
         site != "helby_sw") %>% #remove place holder site for previous figure
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
  mutate(site = case_when(str_detect(transect_id, "EKE02") ~"edward_king_south", #change duplicate name of site that doesnt match up with any other sites
                          TRUE ~ site))%>%
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

#fix helby north lat long error
KCCA_coords <- KCCA_coords %>%
  mutate(latitude = if_else(site == "north_helby_rock", 48.855087, latitude),
         longitude = if_else(site == "north_helby_rock", -125.166609, longitude))

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
         site != "tzartus_116",        #filter these sites out as replaced in 2022
         site != "taylor_island_2021",
         site != "flemming_112_2021") 

# CHECK: NEED TO ASK SANDIE FOR GPS OF NEW SITES AND GENERAL RUN DOWN OF SURVEY AREA ETC

## Combined coords ----
coords <- rbind(JW_coords, DFO_coords, RLS_coords, KRI_coords, SC_coords, KCCA_coords)

### Map ----
# Load GREAT shapefile 
hakai_map <- sf::st_read("data-raw/Hakai_coast/COAST_TEST2.shp") %>%
  st_sf() %>%
  st_set_crs(4326)

# honestly not 100% sure what this does
sf_use_s2(FALSE)

# Define colours
blue <- paste("#b9d1df", sep="")

# transform coords into actually map-able coords
coords_clean <- coords %>%
  st_as_sf(coords = c("longitude", "latitude"),
           remove = FALSE) %>%
  st_set_crs(4326)%>%
  mutate(
    site_label = stringr::str_to_title(gsub("_", " ", site)))

###CHECK WHAT SITES ARE CLOSE TO EACHOTHER AND ON THE SAME ASPECT###


# ============================================================
# IDENTIFY POTENTIALLY MATCHING SITES WITHIN 200 METRES
# ACCOUNTING FOR LAND BARRIERS
# ============================================================

# ------------------------------------------------------------
# 1. PROJECT COORDINATES INTO METRES
# ------------------------------------------------------------

# Project site observations
coords_projected <- coords_clean %>%
  st_transform(32610)

# Project coastline / land polygon
hakai_projected <- hakai_map %>%
  st_transform(32610)

# Combine land polygons into one geometry
land_union <- st_union(hakai_projected)


# ------------------------------------------------------------
# 2. IDENTIFY OBSERVATIONS WITHIN 200 METRES
# ------------------------------------------------------------

nearby_sites <- st_is_within_distance(
  coords_projected,
  coords_projected,
  dist = 230 #This is roughly the size of wizard north and also happens to be the size where sites that are named differently start to be included in the same site df 
)

# Create dataframe of unique site pairs
nearby_pairs <- data.frame(
  pair = rep(
    seq_along(nearby_sites),
    lengths(nearby_sites)
  ),
  match = unlist(nearby_sites)
) %>%
  filter(pair < match) %>%
  mutate(
    site_1 = coords_clean$site[pair],
    site_2 = coords_clean$site[match],
    source_1 = coords_clean$source[pair],
    source_2 = coords_clean$source[match]
  )


# ------------------------------------------------------------
# 3. CALCULATE DISTANCE BETWEEN EACH PAIR
# ------------------------------------------------------------

if (nrow(nearby_pairs) > 0) {
  
  nearby_pairs <- nearby_pairs %>%
    mutate(
      distance_m = as.numeric(
        st_distance(
          coords_projected[pair, ],
          coords_projected[match, ],
          by_element = TRUE
        )
      )
    )
  
} else {
  
  nearby_pairs$distance_m <- numeric(0)
  
}


# ------------------------------------------------------------
# 4. CREATE CONNECTING LINES BETWEEN SITE PAIRS
# ------------------------------------------------------------

if (nrow(nearby_pairs) > 0) {
  
  xy <- st_coordinates(coords_projected)
  
  site_lines <- st_sfc(
    lapply(seq_len(nrow(nearby_pairs)), function(i) {
      
      st_linestring(
        rbind(
          xy[nearby_pairs$pair[i], ],
          xy[nearby_pairs$match[i], ]
        )
      )
      
    }),
    crs = st_crs(coords_projected)
  )
  
  
  # ----------------------------------------------------------
  # 5. CALCULATE HOW MUCH LAND EACH LINE CROSSES
  # ----------------------------------------------------------
  
  land_crossing <- vapply(
    seq_along(site_lines),
    function(i) {
      
      intersection <- st_intersection(
        site_lines[i],
        land_union
      )
      
      if (length(intersection) == 0) {
        return(0)
      }
      
      sum(as.numeric(st_length(intersection)))
      
    },
    numeric(1)
  )
  
  
  # ----------------------------------------------------------
  # 6. CLASSIFY SITE PAIRS
  # ----------------------------------------------------------
  
  nearby_pairs <- nearby_pairs %>%
    mutate(
      land_crossing_m = land_crossing,
      
      # Less than 5 m of land crossing = candidate
      # for being on the same side
      same_side = land_crossing_m < 5
    )
  
} else {
  
  nearby_pairs <- nearby_pairs %>%
    mutate(
      land_crossing_m = numeric(0),
      same_side = logical(0)
    )
  
}


# ------------------------------------------------------------
# 7. SEPARATE POTENTIALLY MATCHING AND SEPARATED SITES
# ------------------------------------------------------------

same_side_sites <- nearby_pairs %>%
  filter(same_side == TRUE) %>%
  arrange(distance_m)

opposite_side_sites <- nearby_pairs %>%
  filter(same_side == FALSE) %>%
  arrange(distance_m)


# ------------------------------------------------------------
# 8. VIEW RESULTS
# ------------------------------------------------------------

# Number of nearby pairs
cat("Total pairs within 200 m:", nrow(nearby_pairs), "\n")

# Number of potentially matching pairs
cat("Potentially matching pairs:", nrow(same_side_sites), "\n")

# Number of pairs separated by land
cat("Pairs crossing land:", nrow(opposite_side_sites), "\n")

# View tables
View(same_side_sites)
View(opposite_side_sites)

##moving same source pairs to opposite_side_sites df
# Identify pairs surveyed by the same source
same_source_pairs <- same_side_sites %>%
  filter(source_1 == source_2)

# Remove these pairs from same_side_sites
same_side_sites <- same_side_sites %>%
  filter(
    source_1 != source_2 |
      grepl("wizard_islet_annual", site_1) |
      grepl("wizard_islet_annual", site_2)
  )

# Add them to opposite_side_sites
opposite_side_sites <- bind_rows(
  opposite_side_sites,
  same_source_pairs
) %>%
  arrange(distance_m)

##manually move sites arranged as different to same that have gps locations too close to land such that the code above assesses them as "oppposite" despite having similar aspect and being of different source origins

# Identify sites of problem to manually return to the "same site" df
manual_assessment <- opposite_side_sites %>%
  filter((grepl("taylor_rock", site_1) | grepl("taylor_rock", site_2)) &
      !(grepl("annual", site_1) | grepl("annual", site_2)) &
        (source_1 != source_2) |
  (grepl("prasiola", site_1) | grepl("prasiola", site_2)) &
      !(grepl("second_beach_south", site_1) | grepl("second_beach_south", site_2)) |
    (grepl("swiss", site_1) | grepl("swiss", site_2)) |
    (grepl("ship", site_1) | grepl("ship", site_2)) |
    (grepl("grappler", site_1) | grepl("grappler", site_2)) |
    (grepl("north_helby_rock", site_1) | grepl("north_helby_rock", site_2)) |
    (grepl("ed_king", site_1) | grepl("ed_king", site_2)) &
    (grepl("edward_king_se", site_1) | grepl("edward_king_se", site_2)))

#add these to same site df
same_side_sites <- bind_rows(
  same_side_sites,
  manual_assessment
) %>%
  distinct(pair, match, .keep_all = TRUE) %>%
  arrange(distance_m)

#remove from opposite site df (if that even matters)
opposite_side_sites <- opposite_side_sites %>%
  anti_join(
    manual_assessment %>%
      dplyr::select(pair, match),
    by = c("pair", "match")
  )

# Create a list of all original observations

# Start with every observation in its own group
coords_grouped <- coords_clean %>%
  mutate(
    observation_id = row_number(),
    location_group = observation_id
  )

# Order approved matching pairs by distance
approved_pairs <- same_side_sites %>%
  arrange(distance_m) %>%
  dplyr::select(pair, match, distance_m)

# Process approved pairs from closest to furthest
for (i in seq_len(nrow(approved_pairs))) {
  
  id_1 <- approved_pairs$pair[i]
  id_2 <- approved_pairs$match[i]
  
  # Identify current groups
  group_1 <- coords_grouped$location_group[id_1]
  group_2 <- coords_grouped$location_group[id_2]
  
  if (group_1 == group_2) next
  
  # Get all observations in each group
  members_1 <- coords_grouped %>%
    filter(location_group == group_1)
  
  members_2 <- coords_grouped %>%
    filter(location_group == group_2)
  
  # Check for duplicate sources
  duplicate_sources <- any(
    members_1$source %in% members_2$source
  )
  
  # Check whether BOTH groups belong to Wizard north
  wizard_exception <- all(
    grepl("wizard",
          c(members_1$site, members_2$site),
          ignore.case = TRUE)
  )
  
  # Merge if sources are unique OR both groups are wizard north
  if (!duplicate_sources || wizard_exception) {
    
    coords_grouped$location_group[
      coords_grouped$location_group == group_2
    ] <- group_1
    
  }
}

# Renumber groups consecutively
coords_grouped <- coords_grouped %>%
  mutate(
    location_group = as.integer(
      factor(location_group)
    )
  )

#check and see if this worked
group_check <- coords_grouped %>%
  st_drop_geometry() %>%
  group_by(location_group) %>%
  summarise(
    n_observations = n(),
    n_sources = n_distinct(source),
    sources = paste(sort(unique(source)), collapse = ", "),
    sites = paste(sort(unique(site)), collapse = ", "),
    .groups = "drop"
  )

View(group_check)

#Nov 6th - make code to join all the names of the grouped sites

# ============================================================
# CREATE ONE CONSOLIDATED COORDINATE PER LOCATION GROUP
# ============================================================

# 1. Project coordinates into metres for distance calculations

coords_projected <- coords_grouped %>%
  st_transform(32610)


# 2. Find the most central observed GPS point in each group
#    (spatial medoid = point with the smallest total distance
#     to all other points in its group)

central_sites <- lapply(
  split(seq_len(nrow(coords_projected)),
        coords_projected$location_group),
  function(i) {
    
    group <- coords_projected[i, ]
    
    distances <- units::drop_units(
      st_distance(group)
    )
    
    group[which.min(rowSums(distances)), ]
    
  }
) %>%
  do.call(what = rbind)


# 3. Keep track of which sources and original site names
#    contributed to each consolidated location

site_metadata <- coords_grouped %>%
  st_drop_geometry() %>%
  group_by(location_group) %>%
  summarise(
    n_observations = n(),
    
    sources = paste(
      sort(unique(source)),
      collapse = ", "
    ),
    
    original_site_names = paste(
      sort(unique(site)),
      collapse = ", "
    ),
    
    .groups = "drop"
  )


# 4. Combine central GPS points with group metadata

coords_consolidated <- central_sites %>%
  rename(
    central_point_source = source
  ) %>%
  left_join(
    site_metadata,
    by = "location_group"
  ) %>%
  st_transform(4326)


# 5. Check results

View(coords_consolidated)

# 6. Order sites by nearest neighbour for mapping 

# ============================================================
# ORDER SITES BY NEAREST NEIGHBOUR
# STARTING AT THE SOUTHERNMOST SITE
# ============================================================

# Project coordinates into metres
coords_ordered <- coords_consolidated %>%
  st_transform(32610)

# Extract coordinates
xy <- st_coordinates(coords_ordered)

# Calculate distance between every pair of sites
dist_matrix <- units::drop_units(
  st_distance(coords_ordered)
)

# Identify southernmost site (smallest Y coordinate)
current <- which.min(xy[, "Y"])

# Create vector to store site numbering
site_order <- integer(nrow(coords_ordered))

# Keep track of sites already visited
visited <- rep(FALSE, nrow(coords_ordered))

# Start numbering at southernmost site
for (i in seq_len(nrow(coords_ordered))) {
  
  site_order[current] <- i
  
  visited[current] <- TRUE
  
  # Stop if all sites have been numbered
  if (i == nrow(coords_ordered)) break
  
  # Find the closest unvisited site
  remaining <- which(!visited)
  
  current <- remaining[
    which.min(dist_matrix[current, remaining])
  ]
}

#check numbering 
site_lookup <- coords_consolidated %>%
  st_drop_geometry() %>%
  dplyr::select(
    site_number,
    location_group,
    site,
    sources,
    longitude,
    latitude
  ) %>%
  arrange(site_number)

View(site_lookup)
# Add numbering back to original dataset
coords_consolidated <- coords_consolidated %>%
  mutate(site_number = site_order)


### END ###

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


## Kieran's minor map revisions
# Legend to the bottom, smaller scale bar, and aligned north arrow
## if you agree, just delete lins 666-762 and it'll flow the same as before

## note- if you end of varying the outplanting site numbers you might need to go with a 3 class system
## note- you've also got some names being double plotting cause they were sampled twice so can switch that easily
# but maybe a version of this map without all site names would be good cause it gets the point across well

## updating
#### Map function ####
map <- function(lat_min, lat_max, long_min, long_max,
                coord_data, nh4_var = NULL, kelp_var, point_size, map_file,
                white_background = TRUE, invert = FALSE) {
  # Invert land and sea colours
  if(invert == FALSE){
    sea <- blue
    land <- "white"
  }
  if(invert == TRUE){
    sea <- "white"
    land <- blue
  }
  # Specify black or white background
  if(white_background == TRUE){
    background <- "white"
    features <- "black"
  }
  if(white_background == FALSE){
    background <- "black"
    features <- "white"
  }
  # Create map
  p <- ggplot() +
    geom_sf(
      data = map_file,
      fill = land,
      colour = sea) +
    # Add points
    # Non-outplant sites (background)
    geom_sf(
      data = coord_data %>% filter(outplant == "no"),
      alpha = 0.9,
      size = point_size,
      stroke = 2,
      aes(
        fill = {{nh4_var}},
        shape = {{kelp_var}},
        colour = {{kelp_var}}
      )
    ) +
    
    # Outplant sites (add second so they are in foreground)
    geom_sf(
      data = coord_data %>% filter(outplant == "yes"),
      alpha = 1,
      size = point_size, #* 2.5, #make outplant sites larger
      stroke = 2,
      aes(
        fill = {{nh4_var}},
        shape = {{kelp_var}},
        colour = {{kelp_var}}
      )
    ) +
    viridis::scale_fill_viridis(
      option = "magma",
      direction = -1,
      guide = guide_colorbar(
        frame.colour = features,
        ticks.colour = features)) +
    coord_sf(
      xlim = c(long_min, long_max),
      ylim = c(lat_min, lat_max),
      expand = FALSE) +
    # ggrepel::geom_text_repel(
    #   data = coord_data,
    #   aes(
    #     x = longitude,
    #     y = latitude
    #     label = site_number),
      # nudge_x = 0.001,
      # direction = "y",
      # hjust = 0,
      # size = 1.5,
      # box.padding = 0.15,
      # point.padding = 0.1,
      # max.overlaps = 5,
      # max.iter = 5000,
      # segment.colour = NA,
      # seed = 123 ) +
      #segment.color = NA) +
    labs(
      fill = "Abalone density",
      shape = "Outplant") +
    scale_shape_manual(
      name = "Outplant",
      values = c(
        "yes" = 21,
        "no" = 21),
      breaks = c("yes", "no")) +
    scale_colour_manual(
      values = c(
        "yes" = "#CD64B5FF",
        "no" = "grey65"),
      breaks = c("yes", "no"),
      guide = "none") +
    guides(
      shape = guide_legend(
        title.position = "left",
        title.hjust = 0.5,
        direction = "horizontal",
        nrow = 1,
        override.aes = list(
          colour = c("#CD64B5FF", "grey65")))) +
    # Themes
    theme_bw() +
    theme(
      panel.background = element_rect(fill = sea),
      panel.grid.major = element_line(color = sea),
      panel.border = element_rect(fill = NA, colour = features),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.ticks.length = grid::unit(0, "pt"),
      plot.margin = grid::unit(c(0, 0, 0, 0), "mm"),
      legend.background = element_rect(color = NA, fill = background),
      legend.key = element_rect(color = background, fill = background),
      legend.key.size = grid::unit(1.2, "lines"),
      legend.text = element_text(size = 24, color = features),
      legend.title = element_text(size = 24, face = "bold", color = features),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(size = 30, color = features),
      strip.text.y = element_text(size = 30, color = features, angle = -90),
      plot.background = element_rect(color = background, fill = background)) +
    # Scale bar - bottom left
    annotation_scale(
      location = "br",
      width_hint = 0.2,
      pad_x = grid::unit(0.1, "in"),
      pad_y = grid::unit(0.05, "in")) +
    # North arrow - top left
    annotation_north_arrow(
      location = "br",
      which_north = "true",
      pad_x = grid::unit(0.01, "in"),
      pad_y = grid::unit(0.15, "in"),
      style = north_arrow_fancy_orienteering)
  return(p)
}

#### Make a map ####
Map_Figure <- map(
  lat_min = 48.73,
  lat_max = 49.00,
  long_min = -125.4,
  long_max = -124.99,
  coord_data = average_abalone_density_site,
  nh4_var = average_density,
  point_size = 3,
  kelp_var = outplant,
  map_file = hakai_map)
Map_Figure

#### Export map at 1200 DPI ####
ggsave(
  filename = "figures/MSc/Abalone_Study_Sites_Map.png",
  plot = Map_Figure,
  width = 20,
  height = 18,
  units = "cm",
  dpi = 1200,
  bg = "white")
##

#rewrite map function with labels for insets
map_label <- function(lat_min, lat_max, long_min, long_max,
                coord_data, nh4_var = NULL, kelp_var, point_size, map_file,
                white_background = TRUE, invert = FALSE) {
  # Invert land and sea colours
  if(invert == FALSE){
    sea <- blue
    land <- "white"
  }
  if(invert == TRUE){
    sea <- "white"
    land <- blue
  }
  # Specify black or white background
  if(white_background == TRUE){
    background <- "white"
    features <- "black"
  }
  if(white_background == FALSE){
    background <- "black"
    features <- "white"
  }
  coord_data <- coord_data %>%
    mutate(
      observation_group = case_when(
        n_observations == 1 ~ "1",
        n_observations == 2 ~ "2",
        n_observations == 3 ~ "3",
        n_observations == 4 ~ "4",
        n_observations == 5 ~ "5",
      )
    )
  # Create map
  p <- ggplot() +
    geom_sf(
      data = map_file,
      fill = land,
      colour = sea) +
    # Add points
    # Non-outplant sites (background)
    geom_sf(
      data = coord_data %>% filter(outplant == "no"),
      alpha = 0.9,
      #size = point_size,
      stroke = 2,
      aes(
        fill = {{nh4_var}},
        shape = observation_group,
        colour = {{kelp_var}},
        size = n_observations
      )
    ) +
    
    # Outplant sites (add second so they are in foreground)
    geom_sf(
      data = coord_data %>% filter(outplant == "yes"),
      alpha = 1,
      #size = point_size * 2.5, #make outplant sites larger
      stroke = 2,
      aes(
        fill = {{nh4_var}},
        shape = observation_group,
        colour = {{kelp_var}},
        size = n_observations
      )
    ) +
    viridis::scale_fill_viridis(
      option = "magma",
      direction = -1,
      guide = guide_colorbar(
        frame.colour = features,
        ticks.colour = features)) +
    scale_size_continuous(
      name = "Number of observations",
      range = c(1, 5)
    ) +
    scale_shape_manual(
      name = "Number of observations",
      values = c(
        "1" = 15,
        "2" = 18,
        "3" = 16,
        "4" = 10,
        "5" = 13
      )
    ) +
    coord_sf(
      xlim = c(long_min, long_max),
      ylim = c(lat_min, lat_max),
      expand = FALSE) +
    # ggrepel::geom_text_repel(
    #   data = coord_data  %>%
    #     filter(
    #       longitude >= long_min,
    #       longitude <= long_max,
    #       latitude >= lat_min,
    #       latitude <= lat_max
    #     ),
    #   aes(
    #     x = longitude,
    #     y = latitude,
    #     label = site_label),
    # nudge_x = 0.001,
    # direction = "y",
    # hjust = 0,
    # size = 1.5,
    # box.padding = 0.15,
    # point.padding = 0.1,
    # max.overlaps = 5,
    # max.iter = 5000,
    # segment.colour = NA,
    # xlim = c(long_min, long_max),
    # ylim = c(lat_min, lat_max),
    # seed = 123 ) +
    labs(
      fill = "Abalone density",
      shape = "Outplant") +
    # #scale_shape_manual(
    #   name = "Outplant",
    #   values = c(
    #     "yes" = 21,
    #     "no" = 21),
    #   breaks = c("yes", "no")) +
    scale_colour_manual(
      values = c(
        "yes" = "#CD64B5FF",
        "no" = "grey65"),
      breaks = c("yes", "no"),
      guide = "none") +
    #guides(
      # shape = guide_legend(
      #   title.position = "left",
      #   title.hjust = 0.5,
      #   direction = "horizontal",
      #   nrow = 1,
      #   override.aes = list(
      #     colour = c("#CD64B5FF", "grey65")))) +
    # Themes
    theme_bw() +
    theme(
      panel.background = element_rect(fill = sea),
      panel.grid.major = element_line(color = sea),
      panel.border = element_rect(fill = NA, colour = features),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.ticks.length = grid::unit(0, "pt"),
      plot.margin = grid::unit(c(0, 0, 0, 0), "mm"),
      legend.background = element_rect(color = NA, fill = background),
      legend.key = element_rect(color = background, fill = background),
      legend.key.size = grid::unit(1.2, "lines"),
      legend.text = element_text(size = 24, color = features),
      legend.title = element_text(size = 24, face = "bold", color = features),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(size = 30, color = features),
      strip.text.y = element_text(size = 30, color = features, angle = -90),
      plot.background = element_rect(color = background, fill = background)) +
    # Scale bar - bottom left
    annotation_scale(
      location = "br",
      width_hint = 0.2,
      pad_x = grid::unit(0.1, "in"),
      pad_y = grid::unit(0.05, "in")) +
    # North arrow - top left
    annotation_north_arrow(
      location = "br",
      which_north = "true",
      pad_x = grid::unit(0.01, "in"),
      pad_y = grid::unit(0.15, "in"),
      style = north_arrow_fancy_orienteering)
  return(p)
}

#SW Deer Group#
SWDEER_Map_Figure <- map_label(
  lat_min = 48.79,
  lat_max = 48.86,
  long_min = -125.26,
  long_max = -125.13,
  coord_data = coords_consolidated,
  # nh4_var = mean_density,
  point_size = 0.5,
  kelp_var = outplant,
  map_file = hakai_map)
SWDEER_Map_Figure

#NE Deer Group#
NEDEER_Map_Figure <- map_label(
  lat_min = 48.8325, #flip these from SW so top = bottom (but add enough to see outplant sites?)
  lat_max = 48.96, 
  long_min = -125.18,
  long_max = -125,
  coord_data = coords_consolidated,
  # nh4_var = mean_density,
  point_size = 0.5,
  kelp_var = outplant,
  map_file = hakai_map)
NEDEER_Map_Figure


## testing the idea of a distance out outplanting grid 
## i'm not sure how you already determined distance to outplanting so don't let this superceed that
## but if you think it's useful then this approach could be revised

library(sf)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggspatial)

#### 1. Prepare spatial data ####
# UTM Zone 10N for distance calculations in metres
sites <- st_transform(coords_clean, 26910)
land <- st_transform(hakai_map, 26910)
outplants <- sites %>% filter(outplant == "yes")

#### 2. Distance from each site to nearest outplant ####
nearest_site <- st_nearest_feature(sites, outplants)
sites$distance_km <- as.numeric(
  st_distance(sites, outplants[nearest_site, ], by_element = TRUE)
) / 1000

#### 3. Create a 250 m grid across the map ####
map_extent <- st_as_sfc(st_bbox(c(
  xmin = -125.4, xmax = -124.99,
  ymin = 48.73, ymax = 49.00
), crs = st_crs(4326))) %>% st_transform(26910)

grid <- st_sf(geometry = st_make_grid(map_extent, cellsize = 100))
grid_centres <- st_centroid(grid)

#### 4. Calculate distance for each grid cell ####
nearest_grid <- st_nearest_feature(grid_centres, outplants)
grid$distance_km <- as.numeric(
  st_distance(grid_centres, outplants[nearest_grid, ], by_element = TRUE)
) / 1000

#### 5. Prepare site labels ####
sites_ll <- st_transform(sites, 4326)
site_labels <- cbind(st_drop_geometry(sites_ll), st_coordinates(sites_ll))

#### 6. Plot distance grid, land, and all sites ####
Distance_Map <- ggplot() +
  # Distance surface
  geom_sf(
    data = grid,
    aes(fill = distance_km),
    colour = NA,
    linewidth = 0) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    name = "Distance to nearest\noutplant (km)") +
  # Land overlay
  geom_sf(
    data = land,
    fill = "white",
    colour = "grey55",
    linewidth = 0.15) +
  # Plot all monitoring sites
  geom_sf(
    data = sites,
    aes(colour = outplant),
    size = 2.5,
    alpha = 0.9) +
  scale_colour_manual(
    values = c("yes" = "#CD64B5FF", "no" = "grey65"),
    name = "Outplant",
    breaks = c("yes", "no")) +
  coord_sf(
    crs = st_crs(4326),
    default_crs = st_crs(4326),
    xlim = c(-125.4, -124.99),
    ylim = c(48.73, 49.00),
    expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 5, "mm"))
Distance_Map

#### 7. Export at 600 DPI ####
ggsave(
  filename = "figures/MSc/Distance_To_Nearest_Outplant.png",
  plot = Distance_Map,
  width = 20,
  height = 18,
  units = "cm",
  dpi = 600,
  bg = "white")

### test end

## AVERAGE DENSITY MAP ----
#FOR MONDAY FIGURE OUT HOW TO COMBINE THE COORDS_CONSOLIDATED DATASET WITH THE AVERAGE DENSITY SCORE BY A SITE, THEN ADD ALL THOSE DENSITIES TOGETHER TO GET AN AVERAGE. MAP WITH EACH POINT COLOUR SCALED BY GRADIENT OF DENSITY TO SEE IF THAT  IS SIMILAR IN A COLOUR PATTERN TO DISTNACE FROM OUTPLANT SITE
###Nov 7th - add an average density score per the 94 sites on the map

# Create lookup table connecting original names to consolidated names
site_lookup <- coords_consolidated %>%
  dplyr::select(site_label, original_site_names, location_group, outplant, geometry) %>%
  distinct() %>%
  separate_rows(original_site_names, sep = ",\\s*") %>%
  mutate(
    original_site_names = str_trim(original_site_names),
    original_site_names = str_to_lower(original_site_names)
  ) %>%
  distinct(original_site_names, .keep_all = TRUE)

View(site_lookup)
#128 observations because some of them directly named the same thing despite different sources

#just use 2021-2026 data (past 5 years)
Average_density_2020s <- read.csv("data-processed/all_abalone_density.csv") %>%
  filter(year > 2020) %>%
  filter(site != "new_taylor_island_2022", #remove for now while dont have GPS!
         site != "new_flemming_112_2022",
         site != "tzartus_116",        #filter these sites out permanently as replaced in 2022 
         site != "taylor_island_2021",
         site != "flemming_112_2021",
         site != "ed_king_?", #site removed for patchy surveying
         site != "self_point") #remove for now while dont have GPS

#join the abalone data, note there should be 345 entries to match all the density scores in that data file
abalone_density_consolidated_coords <- Average_density_2020s %>%
  mutate(
    site_original = site,
    site_match = str_to_lower(str_trim(site))
  ) %>%
  left_join(
    site_lookup,
    by = c("site_match" = "original_site_names")
  ) %>%
  mutate(
    site = coalesce(site_label, site_original)
  ) %>%
  dplyr::select(-site_match, -site_label)%>%
  st_as_sf()

View(abalone_density_consolidated_coords)
#they match 345 observations so all data points are accounted for with GPS locations

#Get an average abalone score across data sets for a given site
# CHECK: THIS IS GOING TO BE MORE HEAVILY WEIGHTED BY RLS DATA AND DFO DATA THAN STORMCOAST AND JW DATA
average_abalone_density_site <- abalone_density_consolidated_coords %>%
  group_by(site, outplant) %>%
  summarise(average_density = mean(site_mean_density), #calculate density average across all years of study across all studies
            number_surveys = n(), #tell us how different surveys go into the average
            number_surveyors = n_distinct(source), #tell us how different surveys go into the average
            .groups = "drop")




  

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





#### Nov 5th: Kieran adding models for analysis of Jane's data ####

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

count(Fig3_glm, outplant,site)  ## confirm that aguilar_point is also now a no

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

## First draft (will replace with model)
# Calculate regression equations for each Site × Period
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

# Give equations different vertical positions
Equation_Data <- Equation_Data %>%
  group_by(site) %>%
  mutate(
    label_y = c(
      max(JaneAnnual$site_mean_density, na.rm = TRUE),
      max(JaneAnnual$site_mean_density, na.rm = TRUE) * 0.90
    )
  ) %>%
  ungroup()

# Plot
LinearRegression_JanesAnnual=ggplot(
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
    x = "Year",
    y =   expression("Density (Abalone per" ~ m^2 * ") <-- Carter confirm this"),
    colour = NULL) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    strip.text = element_text(
      size = 11,
      face = "bold"),
    legend.position = "bottom")

#### Janes Annual Data, Model ####
str(JaneAnnual)

library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(emmeans)

# Keep the two sites and clean factor levels
Model_Data <- JaneAnnual %>%
  filter(site %in% c("taylor", "wizard"), !is.na(site_mean_density), !is.na(year), !is.na(pre_post)) %>%
  droplevels() %>%
  mutate(site = factor(site, levels = c("taylor", "wizard")), pre_post = factor(pre_post, levels = c("pre", "post")))

# Find the first post year separately for each site
Change_Years <- Model_Data %>%
  filter(pre_post == "post") %>%
  group_by(site) %>%
  summarise(change_year = min(year), .groups = "drop")

# Centre year around the pre/post transition for each site
Model_Data <- Model_Data %>%
  left_join(Change_Years, by = "site") %>%
  mutate(year_rel = year - change_year)

# Check numbers of observations
count(Model_Data, site, pre_post)  ## carter confirm the minor difference makes sense

# Fit one model with both intercepts and slopes to differ by site and pre/post
Density_Model <- glmmTMB(site_mean_density ~ site * pre_post * year_rel, family = tweedie(link = "log"), data = Model_Data)
summary(Density_Model)

# Model diagnostics
Density_Model_Residuals <- simulateResiduals(Density_Model, n = 1000)
plot(Density_Model_Residuals) ## looks really good
testDispersion(Density_Model_Residuals)
testZeroInflation(Density_Model_Residuals)
# Compare pre and post slopes within each site
Slope_Comparison <- emtrends(Density_Model, pairwise ~ pre_post | site, var = "year_rel")
Slope_Comparison
Slope_Comparison$contrasts
# Estimated pre/post values at the transition year
Transition_Comparison <- emmeans(Density_Model, pairwise ~ pre_post | site, at = list(year_rel = 0), type = "response")
Transition_Comparison
# Create prediction ranges separately for each Site x Pre/Post combination
Pred_Data <- Model_Data %>%
  group_by(site, pre_post) %>%
  summarise(year_min = min(year), year_max = max(year), .groups = "drop") %>%
  rowwise() %>%
  mutate(year = list(seq(year_min, year_max, length.out = 100))) %>%
  unnest(year) %>%
  dplyr::select(-year_min, -year_max) %>%
  left_join(Change_Years, by = "site") %>%
  mutate(year_rel = year - change_year, site = factor(site, levels = levels(Model_Data$site)), pre_post = factor(pre_post, levels = levels(Model_Data$pre_post)))
# Predict on the model's link scale
Predictions <- predict(Density_Model, newdata = Pred_Data, type = "link", se.fit = TRUE)
# Back-transform predictions and 95% CI to the density scale
Pred_Data <- Pred_Data %>%
  mutate(pred = exp(Predictions$fit), lower = exp(Predictions$fit - 1.96 * Predictions$se.fit), upper = exp(Predictions$fit + 1.96 * Predictions$se.fit))
# View the actual predicted dataset
Pred_Data

## We are going to plot some information on the graph
# currently it's the SARA line, but anything can be added the same way

# Plot raw data + model predictions + 95% CI
SARA_Line1 <- data.frame(site = unique(Model_Data$site), x = 2003, xend = 2003, y = 0, yend = 1.2, label_y = 1.2)
Density_Figure <- ggplot() +
  geom_ribbon(data = Pred_Data, aes(x = year, ymin = lower, ymax = upper, fill = pre_post), alpha = 0.20, colour = NA) +
  geom_line(data = Pred_Data, aes(x = year, y = pred, colour = pre_post), linewidth = 1.1) +
  geom_point(data = Model_Data, aes(x = year, y = site_mean_density, colour = pre_post), size = 2.5, alpha = 0.75) +
  geom_segment(data = SARA_Line1, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, linewidth = 0.6, colour = "black") +
  geom_text(data = SARA_Line1, aes(x = x, y = label_y, label = "SARA Listed"), inherit.aes = FALSE, angle = 0, vjust = -0.15, hjust = 0.5, size = 3.5) +
  facet_wrap(~site, nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_colour_manual(values = c("pre" = "#F8766D", "post" = "#00BFC4"), labels = c("pre" = "Pre-outplanting", "post" = "Post-outplanting"), name = NULL) +
  scale_fill_manual(values = c("pre" = "#F8766D", "post" = "#00BFC4"), labels = c("pre" = "Pre-outplanting", "post" = "Post-outplanting"), name = NULL) +
  labs(x = "Year", y = expression("Density (Abalone per" ~ m^2 * ")"), colour = NULL, fill = NULL) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), strip.text = element_text(size = 11, face = "bold"), legend.position = "bottom")
Density_Figure


#### If Isa wants predicted data points 
# Generate one fitted/predicted value for every raw observation
Model_Data <- Model_Data %>%
  mutate(predicted_density = predict(Density_Model, newdata = Model_Data, type = "response"),
         observation = row_number())
# Check that predicted and observed have the same number of values
nrow(Model_Data)
length(Model_Data$predicted_density)
# Compare observed and predicted values
t.test(Model_Data$site_mean_density, Model_Data$predicted_density, paired = TRUE)
## p-value is 0.97 so the take home is no difference between raw and predicted

# One important note, technically, the predictions above a fitted expected values (so like given model predict this)
# We can also generate a synthetic dataset 
set.seed(123)
Sim_Data <- simulate(Density_Model, nsim = 1)
Model_Data$simulated_density <- Sim_Data$sim_1
Model_Data$site_mean_density
t.test(Model_Data$site_mean_density, Model_Data$simulated_density, paired = TRUE) ## again we'll test it
## we'll also plot the t.test data
## p-value is 0.94. 
DataPredict <- Model_Data %>%
  mutate(observation = row_number()) %>%
  dplyr::select(observation, site, pre_post, site_mean_density, simulated_density) %>%
  pivot_longer(cols = c(site_mean_density, simulated_density), names_to = "type", values_to = "density") %>%
  mutate(type = recode(type, site_mean_density = "Observed", simulated_density = "Simulated"))
Observed_Simulated_Plot <- ggplot(DataPredict, aes(x = observation, y = density, colour = type)) +
  geom_point(size = 2, shape = 1) +
  scale_colour_manual(values = c("Observed" = "red", "Simulated" = "black")) +
  theme_classic() +
  theme(legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 13, colour = "black"), axis.text = element_text(colour = "black"), text = element_text(size = 15, colour = "black")) +
  labs(x = "Observation", y = expression("Density (Abalone per" ~ m^2 * ")"))
# I feel like a take home of this approach is that raw date, predictions, expected, etc. should be similar if you do the modelling well


## we can then plot these
SARA_Line <- data.frame(site = c("taylor", "wizard"), x = c(2003, 2003), xend = c(2003, 2003), y = c(0, 0), yend = c(1.2, 1.2), label_y = c(2.08, 2.08))

Simulated_Figure <- ggplot() +
  geom_ribbon(data = Pred_Data, aes(x = year, ymin = lower, ymax = upper, fill = pre_post), alpha = 0.20, colour = NA) +
  geom_line(data = Pred_Data, aes(x = year, y = pred, colour = pre_post), linewidth = 1.1) +
  geom_point(data = Model_Data, aes(x = year, y = simulated_density, colour = pre_post), size = 2.5, alpha = 0.8) +
  geom_segment(data = SARA_Line, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, linewidth = 0.6, colour = "black") +
  geom_text(data = SARA_Line, aes(x = x, y = 1.2, label = "SARA Listed"), inherit.aes = FALSE, angle = 0, vjust = -0.2, hjust = 0.5, size = 3.5) +
  facet_wrap(~site, nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_colour_manual(values = c("pre" = "#F8766D", "post" = "#00BFC4"), labels = c("pre" = "Pre-outplanting", "post" = "Post-outplanting")) +
  scale_fill_manual(values = c("pre" = "#F8766D", "post" = "#00BFC4"), guide = "none") +
  labs(x = "Year", y = expression("Density (Abalone per" ~ m^2 * ")"), colour = NULL) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), strip.text = element_text(size = 11, face = "bold"), legend.position = "bottom")
Simulated_Figure


## Summary
LinearRegression_JanesAnnual ## initial linear model
Density_Figure   ## modeled values with raw data plotted
Simulated_Figure ## modeled values with prediced data plotted

## The model is valid 
plot(Density_Model_Residuals)
testDispersion(Density_Model_Residuals)
testZeroInflation(Density_Model_Residuals)
t.test(Model_Data$site_mean_density, Model_Data$simulated_density, paired = TRUE) ## again we'll test it
Observed_Simulated_Plot

plot(Density_Model_Residuals)
testDispersion(Density_Model_Residuals)
testZeroInflation(Density_Model_Residuals)
Observed_Simulated_Plot

# Exporting
## if you use this, you'll want to tinker the dimensions a bit cause I'm just guessing
Diagnostics_Figure_2x2 <- DHARMa_Plot + Dispersion_Plot + ZeroInflation_Plot + Observed_Simulated_Plot
ggsave("figures/MSc/Model_Diagnostics_2x2.png", Diagnostics_Figure_2x2, width = 18, height = 18, units = "in", dpi = 600)

Simulated_Figure
ggsave("figures/MSc/JaneAnnualSurvey_Simulated_Figure.png", Simulated_Figure, width = 9, height = 5, units = "in", dpi = 600)







#### Jane Random Surveys #### 
JaneRandom=Fig3_glm[Fig3_glm$method == "par_quad", ]
count(JaneRandom,site,outplant)
str(JaneRandom)
hist(JaneRandom$site_mean_density)

### Inital examination Nov 1st 2026
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

#### Nov 5 Kieran Addition, model of Jane Random Data ####
str(JaneRandom)
outplant
count(JaneRandom,year,outplant)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(dplyr)
library(tidyr)
library(ggplot2)
JaneRandom <- JaneRandom %>%
  droplevels() %>%
  mutate(year = factor(year, levels = c("1988", "1994", "2007", "2022")), outplant = factor(outplant, levels = c("no", "yes")))
Density_Model_Jane <- glmmTMB(site_mean_density ~ year * outplant + (1 | site), family = tweedie(link = "log"), data = JaneRandom)
summary(Density_Model_Jane)

## testing the quality of the model
Jane_Residuals <- simulateResiduals(Density_Model_Jane, n = 1000)
plot(Jane_Residuals)
testDispersion(Jane_Residuals)
testZeroInflation(Jane_Residuals)

Jane_Residuals <- simulateResiduals(Density_Model_Jane, n = 1000)
plot(Jane_Residuals)
testDispersion(Jane_Residuals)
testZeroInflation(Jane_Residuals)
# Generate one synthetic dataset from the fitted model
set.seed(123)
Jane_Sim <- simulate(Density_Model_Jane, nsim = 1)
JaneRandom$simulated_density <- Jane_Sim$sim_1
# Check that observed and simulated datasets have the same number of values
nrow(JaneRandom)
length(JaneRandom$simulated_density)
# Compare observed and simulated values
Jane_Sim_ttest <- t.test(JaneRandom$site_mean_density, JaneRandom$simulated_density, paired = TRUE)
Jane_Sim_ttest ## again, no sign difference between raw and predicted data

## let's plot it
Jane_DataPredict <- JaneRandom %>%
  mutate(observation = row_number()) %>%
  dplyr::select(observation, year, outplant, site_mean_density, simulated_density) %>%
  pivot_longer(cols = c(site_mean_density, simulated_density), names_to = "type", values_to = "density") %>%
  mutate(type = recode(type, site_mean_density = "Observed", simulated_density = "Simulated"))
Jane_Observed_Simulated_Plot <- ggplot(Jane_DataPredict, aes(x = observation, y = density, colour = type)) +
  geom_point(size = 2, shape = 1) +
  scale_colour_manual(values = c("Observed" = "red", "Simulated" = "black")) +
  theme_classic() +
  theme(legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 13, colour = "black"), axis.text = element_text(colour = "black"), text = element_text(size = 15, colour = "black")) +
  labs(x = "Observation", y = expression("Abalone density (" * m^-2 * ")"))
Jane_Observed_Simulated_Plot

## running stats on this data
Jane_EMM <- emmeans(Density_Model_Jane, ~ outplant | year, type = "response")
Jane_Pred <- as.data.frame(Jane_EMM)
Jane_Pred
Jane_Contrasts <- pairs(Jane_EMM, adjust = "holm")
Jane_Contrasts ## 2022 is sign different, while all other years aren't

Jane_EMM_Year <- emmeans(Density_Model_Jane, ~ year | outplant, type = "response")
Jane_Year_Contrasts <- pairs(Jane_EMM_Year, adjust = "tukey")
Jane_Year_Contrasts ## this looks at within outplant status 

Jane_Pred_Figure <- ggplot() +
  geom_point(data = JaneRandom, aes(x = year, y = site_mean_density, colour = outplant), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.35), size = 2, alpha = 0.45) +
  geom_errorbar(data = Jane_Pred, aes(x = year, ymin = asymp.LCL, ymax = asymp.UCL, colour = outplant, group = outplant), position = position_dodge(width = 0.35), width = 0.12, linewidth = 0.7) +
  geom_point(data = Jane_Pred, aes(x = year, y = response, colour = outplant, group = outplant), position = position_dodge(width = 0.35), size = 3) +
  scale_colour_manual(values = c("no" = "#F8766D", "yes" = "#00BFC4"), labels = c("no" = "No", "yes" = "Yes")) +
  labs(x = "Year", y = expression("Abalone density (" * m^-2 * ")"), colour = "Outplant status") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), legend.position = "bottom")
Jane_Pred_Figure

## matrix of linear relationships ##
str(JaneRandom)

## we'll do this as a first pass then come back to check the model in detail
# I'm worried about so many predictions with min data (but lets see how it fits first)

# numeric year and centre it for the model
JaneRandom <- JaneRandom %>%
  mutate(year_num = as.numeric(as.character(year)),
         year_centered = (year_num - 2005) / 10)
# Fit model: overall effect of outplant status, with site-specific intercepts and slopes
Jane_Site_Model <- glmmTMB(site_mean_density ~ year_centered * outplant + (year_centered | site),
                           family = tweedie(link = "log"),
                           data = JaneRandom)
summary(Jane_Site_Model)
# Diagnostics
Jane_Site_Residuals <- simulateResiduals(Jane_Site_Model, n = 1000)
plot(Jane_Site_Residuals)  ## much better then I thought it would be
testDispersion(Jane_Site_Residuals)
testZeroInflation(Jane_Site_Residuals)
# Create a continuous prediction sequence for every site
Jane_Site_Pred <- JaneRandom %>%
  dplyr::select(site, outplant) %>%
  distinct() %>%
  tidyr::crossing(year_num = seq(1988, 2022, length.out = 100)) %>%
  mutate(year_centered = (year_num - 2005) / 10)
# Generate site-specific predictions INCLUDING random effects
Jane_Site_Predict <- predict(Jane_Site_Model,
                             newdata = Jane_Site_Pred,
                             type = "link",
                             se.fit = TRUE,
                             re.form = NULL)
# Back-transform predictions and confidence intervals
Jane_Site_Pred <- Jane_Site_Pred %>%
  mutate(pred = exp(Jane_Site_Predict$fit),
         lower = exp(Jane_Site_Predict$fit - 1.96 * Jane_Site_Predict$se.fit),
         upper = exp(Jane_Site_Predict$fit + 1.96 * Jane_Site_Predict$se.fit))
# Plot site-specific predictions + raw observations
Jane_Site_Figure <- ggplot() +
  geom_ribbon(data = Jane_Site_Pred, aes(x = year_num, ymin = lower, ymax = upper, fill = outplant), alpha = 0.20, colour = NA) +
  geom_line(data = Jane_Site_Pred, aes(x = year_num, y = pred, colour = outplant), linewidth = 1.1) +
  geom_point(data = JaneRandom, aes(x = year_num, y = site_mean_density, colour = outplant), size = 2.5, alpha = 0.8) +
  facet_wrap(~site, ncol = 5) +
  scale_colour_manual(values = c("yes" = "#CD64B5", "no" = "grey60"),
                      labels = c("yes" = "yes", "no" = "no")) +
  scale_fill_manual(values = c("yes" = "#CD64B5", "no" = "grey60"),
                    guide = "none") +
  scale_x_continuous(breaks = c(1988, 2007, 2022)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "", y = expression("Abalone density (" * m^-2 * ")"), colour = "Outplant") +
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 9),
        legend.position = "right")
Jane_Site_Figure
##

ggsave(
  filename = "figures/MSc/Abalone_Jane_Site_Figure.png",
  plot = Jane_Site_Figure,
  width = 20,
  height = 16,
  units = "cm",
  dpi = 900,
  bg = "white")

## note, the points are raw not predicted data but predictions are easy to add
# although I think in this case (and really most cases) it's nice to have raw data
# cause you then see if your model matches your raw data, if you have predictions 
# you just see if your predictions match your prediction ahh





##






#







#

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

