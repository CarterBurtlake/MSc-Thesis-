# 0a. Header & Metadata ----------------------------------------------------
# Title: [MSc / Figures / Clean]
# Author: [Carter Burtlake]
# Date: [2026-11-24]
# Purpose: [generates the figures for thesis chapter 2]

# 0b. Setup & Libraries ----------------------------------------------------
library(readr)
library(tidyr)
library(purrr) # for functions and loops (bootstrapping)
library(dplyr) # for general data wrangling
library(lubridate) # for looking at dates
library(janitor) #to clean data frames
library(stringr) #to clean data frames
library(forcats) #to order factors
library(tibble) #for tabulated data

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

# 1. Fig 1: Map -----------------------------------------------

# 2. Fig 2: Graphical representation of map -----------------------------------------------

# 3. Fig 3: JW data isolated -----------------------------------------------

#first make df have site level outplant information
rate <- read.csv("data-processed/all_abalone_density.csv") %>%
  mutate(outplant = case_when(str_detect(site, "scotts_bay") ~"yes", #add a column for outplant level
                              str_detect(site, "aguilar_point") ~"yes",
                              str_detect(site, "grappler") ~"yes",
                              str_detect(site, "taylor_ed_king") ~"yes",
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
                         "aguilar_point",
                         "execution_rock",
                         "wizard",
                         "helby_island",
                         "self_point",
                         "whittlestone",
                         "prasiola",
                         "cape_beale",
                         "w_ed_king",
                         "ship_islands",
                         "seppings",
                         "lawton_point",
                         "cia_rock",
                         "kirby_point",
                         "village_bay",
                         "blowhole",
                         "taylor_ed_king",
                         "ed_king",
                         "taylor")))

Fig3_glm <- Fig3_rate %>%
  #filter(method != "perp_belt") %>% #remove annual data for modelling 
  #notice that if you add the filter above you take a lot of the weight out of the non_outplant data that is rich in the annual surveys and this drives a non-significant pattern with the interaction but generates one with year
  mutate(year_decade = (year - mean(year, na.rm = TRUE)) / 10) %>% #attempt to make year continuous for the random sites and see "change per decade or per year"
  mutate(pre_post = case_when(year > 2003 ~ "post",
                              TRUE ~ "pre")) %>%                    #set a pre outplant vs post outplant time period
  mutate(outplant = case_when(site == "taylor_ed_king" ~ "no",
                              site == "aguilar_point" ~ "no",
                              TRUE ~ outplant))

#annual sites
JaneAnnual=Fig3_glm[Fig3_glm$method == "perp_belt", ]

#add factor and make sure Pre/Post is ordered correctly
levels(as.factor(JaneAnnual$pre_post))
JaneAnnual$pre_post=as.factor(JaneAnnual$pre_post)
JaneAnnual$pre_post <- factor(
  JaneAnnual$pre_post,
  levels = c("pre", "post"))

####assess linear relationships for the annual sites####
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

LinearRegression_JanesAnnual
#CHECK: why is the R squared value for pre wizard 0...

#####Model annual sites####
# Keep the two sites and clean factor levels
Model_Data_JW_Annual <- JaneAnnual %>%
  filter(site %in% c("taylor", "wizard"), !is.na(site_mean_density), !is.na(year), !is.na(pre_post)) %>%
  droplevels() %>%
  mutate(site = factor(site, levels = c("taylor", "wizard")), pre_post = factor(pre_post, levels = c("pre", "post")))

# Find the first post year separately for each site
Change_Years <- Model_Data_JW_Annual %>%
  filter(pre_post == "post") %>%
  group_by(site) %>%
  summarise(change_year = min(year), .groups = "drop")

# Center year around the pre/post transition for each site
Model_Data_JW_Annual <- Model_Data_JW_Annual %>%
  left_join(Change_Years, by = "site") %>%
  mutate(year_rel = year - change_year)

# Check numbers of observations
count(Model_Data_JW_Annual, site, pre_post)  ## Taylor missing site recording in 1991 as expected per comms JW

# Fit one model with both intercepts and slopes to differ by site and the relative year to outplanting
Annual_Density_Model <- glmmTMB(site_mean_density ~ site * year_rel *pre_post, family = tweedie(link = "log"), data = Model_Data_JW_Annual)
###CHECK: GET KIERAN TO VERIFY THE INCLUSION OF BOTH TERMS (YEAR REL IS INCLUDED AS IT ENSURES THE MODEL IS INTERPRETTING THE YEARS CORRECTLY OVER THE NUMBER LINE, AND PRE POST LETS US ASSESS WHETHER IT THERE IS AN OUTPLANT EFFECT). Maybe move site to RE?
summary(Annual_Density_Model)

# Model diagnostics
Annual_Density_Model_Residuals <- simulateResiduals(Annual_Density_Model, n = 1000)
plot(Annual_Density_Model_Residuals) ## looks really good
testDispersion(Annual_Density_Model_Residuals)
testZeroInflation(Annual_Density_Model_Residuals)
# Compare pre and post slopes within each site
Annual_JW_Slope_Comparison <- emtrends(Annual_Density_Model, pairwise ~ pre_post | site, var = "year_rel")
Annual_JW_Slope_Comparison
Annual_JW_Slope_Comparison$contrasts
# Estimated pre/post values at the transition year
Annual_JW_Transition_Comparison <- emmeans(Annual_Density_Model, pairwise ~ pre_post | site, at = list(year_rel = 0), type = "response")
Annual_JW_Transition_Comparison
# Create prediction ranges separately for each Site x Pre/Post combination
Annual_JW_Pred_Data <- Model_Data_JW_Annual %>%
  group_by(site, pre_post) %>%
  summarise(year_min = min(year), year_max = max(year), .groups = "drop") %>%
  rowwise() %>%
  mutate(year = list(seq(year_min, year_max, length.out = 100))) %>%
  unnest(year) %>%
  dplyr::select(-year_min, -year_max) %>%
  left_join(Change_Years, by = "site") %>%
  mutate(year_rel = year - change_year, site = factor(site, levels = levels(Model_Data_JW_Annual$site)), pre_post = factor(pre_post, levels = levels(Model_Data_JW_Annual$pre_post)))
# Predict on the model's link scale
Annual_JW_Predictions <- predict(Annual_Density_Model, newdata = Annual_JW_Pred_Data, type = "link", se.fit = TRUE)
# Back-transform predictions and 95% CI to the density scale
Annual_JW_Pred_Data <- Annual_JW_Pred_Data %>%
  mutate(pred = exp(Annual_JW_Predictions$fit), lower = exp(Annual_JW_Predictions$fit - 1.96 * Annual_JW_Predictions$se.fit), upper = exp(Annual_JW_Predictions$fit + 1.96 * Annual_JW_Predictions$se.fit))
# View the actual predicted dataset
Annual_JW_Pred_Data

# Plot raw data + model predictions + 95% CI
SARA_Line <- data.frame(site = unique(Model_Data_JW_Annual$site), x = 2003, xend = 2003, y = 0, yend = 1.2, label_y = 1.2)
#add when the species was SARA listed as threatened
Harvest_Line <- data.frame(site = unique(Model_Data_JW_Annual$site), x = 1990, xend = 1990, y = 0, yend = 1.2, label_y = 1.2)
#add when the species was prohibitted from harvest
Annual_JW_Density_Figure <- ggplot() +
  geom_ribbon(data = Annual_JW_Pred_Data, aes(x = year, ymin = lower, ymax = upper, fill = pre_post), alpha = 0.20, colour = NA) +
  geom_line(data = Annual_JW_Pred_Data, aes(x = year, y = pred, colour = pre_post), linewidth = 1.1) +
  geom_point(data = Model_Data_JW_Annual, aes(x = year, y = site_mean_density, colour = pre_post), size = 2.5, alpha = 0.75) +
  annotate("rect",
           xmin = 2003, xmax = 2010,
           ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.3) +
  annotate("text",
           x = 2006.5, y = 2.1,
           label = "Outplanting",
           vjust = 1.5, size = 2.5)+
  geom_segment(data = SARA_Line, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, linewidth = 0.6, colour = "black") +
  geom_segment(data = Harvest_Line, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, linewidth = 0.6, colour = "black") +
  geom_text(data = SARA_Line, aes(x = x, y = label_y, label = "SARA listed"), inherit.aes = FALSE, angle = 0, vjust = -0.15, hjust = 0.5, size = 2.5) +
  geom_text(data = Harvest_Line, aes(x = x, y = label_y, label = "Harvest ban"), inherit.aes = FALSE, angle = 0, vjust = -0.15, hjust = 0.5, size = 2.5) +
  facet_wrap(~site, nrow = 1,
             labeller = as_labeller(c("taylor" = "Taylor",
                                      "wizard" = "Wizard"))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_colour_manual(values = c("pre" = "#F8766D", "post" = "#00BFC4"), labels = c("pre" = "Pre-outplanting", "post" = "Post-outplanting"), name = NULL) +
  scale_fill_manual(values = c("pre" = "#F8766D", "post" = "#00BFC4"), labels = c("pre" = "Pre-outplanting", "post" = "Post-outplanting"), name = NULL) +
  labs(x = "Year", y = expression("Density (Abalone per" ~ m^2 * ")"), colour = NULL, fill = NULL) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), strip.text = element_text(size = 11, face = "bold"), legend.position = "bottom")
Annual_JW_Density_Figure

#can include simulated data based of the model, but I don't think that is necessary

## Summary
LinearRegression_JanesAnnual ## initial linear model
Annual_JW_Density_Figure ## raw data with CI calculated using predict function

##JW Random data ----------
#main outplant sites df
JaneRandom=Fig3_glm[Fig3_glm$method == "par_quad", ]
hist(JaneRandom$site_mean_density)
JaneRandom$year=as.factor(JaneRandom$year)

#assign outplant levels as factor for graphing and more useful model comparisons
#NOTE IF YOU SHIFT THIS THE MODEL OUTPUT FOR YEAR CHANGES AS THE INTERCEPT EITHER COMPARES NON-OUTPLANTS ACROSS YEARS OR OUTPLANTS ACROSS YEARS
#if you dont run this code the model reports a significant -ve decline in density at non-outplant sites in 2007 relative to 1988. If you do it reports a significant +ve increase in density at outplant sites in 2022 relative to 1988 (despite the average effect being non-significant / no change across all years)
JaneRandom<- JaneRandom %>%
  mutate(outplant = factor(outplant, levels = c("yes", "no")))

#make df without annual but with medium outplant sites
medium_outplant_JW <- Fig3_rate %>%
  #filter(method != "perp_belt") %>% #remove annual data for modelling 
  #notice that if you add the filter above you take a lot of the weight out of the non_outplant data that is rich in the annual surveys and this drives a non-significant pattern with the interaction but generates one with year
  mutate(year_decade = (year - mean(year, na.rm = TRUE)) / 10) %>% #attempt to make year continuous and see "change per decade"
  mutate(pre_post = case_when(year > 2003 ~ "post",
                              TRUE ~ "pre")) %>%
  filter(method == "par_quad")

#assign year as factor 
medium_outplant_JW$year=as.factor(medium_outplant_JW$year)

#assign outplant levels as factor for graphing
medium_outplant_JW<- medium_outplant_JW %>%
  mutate(outplant = factor(outplant, levels = c("yes", "no")))

####model year as factor based off main outplants####
Main_JaneRandom_Model <- glmmTMB(site_mean_density ~ as.factor(year)*outplant+ (1| site), family = tweedie(link = "log"), data = JaneRandom)
plot(simulateResiduals(Main_JaneRandom_Model)) #looks good
summary(Main_JaneRandom_Model) #looks to be an effect of recent years with +ve abalone density

## testing the quality of the model
Main_Jane_Residuals <- simulateResiduals(Main_JaneRandom_Model, n = 1000)
plot(Main_Jane_Residuals)
testDispersion(Main_Jane_Residuals)
testZeroInflation(Main_Jane_Residuals)

## lets generate a plot to visually see this effect
# generate synthetic dataset from fitted model to get summary stats (and report error accurately)
set.seed(123)
Jane_Sim1 <- simulate(Main_JaneRandom_Model, nsim = 1)
JaneRandom$simulated_density <- Jane_Sim1$sim_1

# Check that observed and simulated data sets have the same number of values
nrow(JaneRandom)
length(JaneRandom$simulated_density) #both 80

# Compare observed and simulated values
Jane_Sim_ttest <- t.test(JaneRandom$site_mean_density, JaneRandom$simulated_density, paired = TRUE)
Jane_Sim_ttest ## no sig difference between raw and predicted data

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
Main_JaneRanom_Emmeans <- emmeans(Main_JaneRandom_Model, ~ outplant | year, type = "response")
Jane_Pred1 <- as.data.frame(Main_JaneRanom_Emmeans)
Jane_Pred1
Jane_Contrasts1 <- pairs(Main_JaneRanom_Emmeans, adjust = "holm")
Jane_Contrasts1 ## 2022 is sign different, while all other years aren't

Main_JaneRanom_Emmeans_year <- emmeans(Main_JaneRandom_Model, ~ year | outplant, type = "response")
Main_Jane_Year_Contrasts <- pairs(Main_JaneRanom_Emmeans_year, adjust = "tukey")
Main_Jane_Year_Contrasts ## this looks at within outplant status 

Main_Jane_Pred_Figure <- ggplot() +
  geom_point(data = JaneRandom, aes(x = year, y = site_mean_density, colour = outplant), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.35, reverse =TRUE), size = 2, alpha = 0.45) +
  geom_errorbar(data = Jane_Pred1, aes(x = year, ymin = asymp.LCL, ymax = asymp.UCL, colour = outplant, group = outplant), position = position_dodge(width = 0.35, reverse = FALSE), width = 0.12, linewidth = 0.7) +
  geom_point(data = Jane_Pred1, aes(x = year, y = response, colour = outplant, group = outplant), position = position_dodge(width = 0.35, reverse = FALSE), size = 3) +
  scale_colour_manual(values = c("no" = "grey60", "yes" = "#CD64B5"), labels = c("no" = "Non-outplant", "yes" = "Outplant")) +
  #guides(color = guide_legend(reverse = TRUE)) + #reverse legend to get yes on top
  labs(x = "Year", y = expression("Abalone density (" * m^-2 * ")"), colour = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), legend.position = "bottom")
Main_Jane_Pred_Figure
#seems like we cant lump pre/post into a category as the 2007 data would average the 2022 data. 
#lets look at this with all ("medium") outplanting sites

####model year as factor based off medium outplants####
Medium_JaneRandom_Model <- glmmTMB(site_mean_density ~ as.factor(year)*outplant+ (1| site), family = tweedie(link = "log"), data = medium_outplant_JW)
plot(simulateResiduals(Medium_JaneRandom_Model)) #looks good
summary(Medium_JaneRandom_Model) #nothing of interest

## testing the quality of the model
Medium_Jane_Residuals <- simulateResiduals(Medium_JaneRandom_Model, n = 1000)
plot(Medium_Jane_Residuals)
testDispersion(Medium_Jane_Residuals)
testZeroInflation(Medium_Jane_Residuals)

## lets generate a plot to visually see this effect
# generate synthetic dataset from fitted model to get summary stats (and report error accurately)
set.seed(123)
Jane_Sim2 <- simulate(Medium_JaneRandom_Model, nsim = 1)
medium_outplant_JW$simulated_density <- Jane_Sim2$sim_1

# Check that observed and simulated data sets have the same number of values
nrow(medium_outplant_JW)
length(medium_outplant_JW$simulated_density) #both 80

# Compare observed and simulated values
Jane_Sim_ttest2 <- t.test(medium_outplant_JW$site_mean_density, medium_outplant_JW$simulated_density, paired = TRUE)
Jane_Sim_ttest2 ## no sig difference between raw and predicted data

## let's plot it
Jane_DataPredict2 <- medium_outplant_JW %>%
  mutate(observation = row_number()) %>%
  dplyr::select(observation, year, outplant, site_mean_density, simulated_density) %>%
  pivot_longer(cols = c(site_mean_density, simulated_density), names_to = "type", values_to = "density") %>%
  mutate(type = recode(type, site_mean_density = "Observed", simulated_density = "Simulated"))

Jane_Observed_Simulated_Plot2 <- ggplot(Jane_DataPredict2, aes(x = observation, y = density, colour = type)) +
  geom_point(size = 2, shape = 1) +
  scale_colour_manual(values = c("Observed" = "red", "Simulated" = "black")) +
  theme_classic() +
  theme(legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 13, colour = "black"), axis.text = element_text(colour = "black"), text = element_text(size = 15, colour = "black")) +
  labs(x = "Observation", y = expression("Abalone density (" * m^-2 * ")"))
Jane_Observed_Simulated_Plot

## running stats on this data
Medium_JaneRanom_Emmeans <- emmeans(Medium_JaneRandom_Model, ~ outplant | year, type = "response")
Jane_Pred2 <- as.data.frame(Medium_JaneRanom_Emmeans)%>%
  mutate(outplant = factor(outplant, levels = c("yes", "no")))
Jane_Pred2
Jane_Contrasts2 <- pairs(Medium_JaneRanom_Emmeans, adjust = "holm")
Jane_Contrasts2 ## 2022 is sign different, while all other years aren't. Maybe moderate evidence for a -ve effect in 1994 at outplant = yes sites

Medium_JaneRanom_Emmeans_year <- emmeans(Medium_JaneRandom_Model, ~ year | outplant, type = "response")
Medium_Jane_Year_Contrasts <- pairs(Medium_JaneRanom_Emmeans_year, adjust = "tukey")
Medium_Jane_Year_Contrasts ## this looks at within outplant status 

Medium_Jane_Pred_Figure <- ggplot() +
  geom_point(data = medium_outplant_JW, aes(x = year, y = site_mean_density, colour = outplant), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.35, reverse = TRUE), size = 2, alpha = 0.45) +
  geom_errorbar(data = Jane_Pred2, aes(x = year, ymin = asymp.LCL, ymax = asymp.UCL, colour = outplant, group = outplant), position = position_dodge(width = 0.35, reverse = FALSE), width = 0.12, linewidth = 0.7) +
  geom_point(data = Jane_Pred2, aes(x = year, y = response, colour = outplant, group = outplant), position = position_dodge(width = 0.35, reverse = FALSE), size = 3) +
  scale_colour_manual(values = c("no" = "grey60", "yes" = "#CD64B5"), labels = c("no" = "Non-outplant", "yes" = "Outplant")) +
  #guides(color = guide_legend(reverse = TRUE)) + #reverse legend to get yes on top
  labs(x = "Year", y = expression("Abalone density (" * m^-2 * ")"), colour = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), legend.position = "bottom")
Medium_Jane_Pred_Figure

####Comparing models to consider rate of recovery####
#IF I WANT TO COMPARE ACROSS SITES I NEED TO HAVE A NUMERICAL YEAR FACTOR THAT IS CENTERED AND CONTINUOUS - Also, to review the rates across sites we need to add a random slope for each site across all x data points i.e., including site and year_centered in the RE
# numeric year and centre it for the model with main outplant sites
JaneRandomYearCentered <- JaneRandom %>%
  mutate(year_num = as.numeric(as.character(year)),
         year_centered = (year_num - min(year_num, na.rm = TRUE)) / 10) #sets intercept as time 1988 (previously had this centered and it showed a nearly significant effect in 2005 for the outplantyes intercept, which would be misleading)
# Fit model: overall effect of outplant status, with site-specific intercepts and slopes
Jane_Site_Model <- glmmTMB(site_mean_density ~ year_centered * outplant + (year_centered | site),
                           family = tweedie(link = "log"),
                           data = JaneRandomYearCentered)
summary(Jane_Site_Model) #output slightly varies with or without year_centered considered as RE, but story of signifigance remains the same

## testing the quality of the model
Jane_Site_Residuals <- simulateResiduals(Jane_Site_Model, n = 1000)
plot(Jane_Site_Residuals) #overall looks ok - good (but not perfect)
testDispersion(Jane_Site_Residuals)
testZeroInflation(Jane_Site_Residuals)

# Fit model: overall effect of outplant status, withOUT site-specific intercepts and slopes
Jane_Site_Model2 <- glmmTMB(site_mean_density ~ year_centered * outplant + (1| site),
                           family = tweedie(link = "log"),
                           data = JaneRandomYearCentered)
summary(Jane_Site_Model2) #output slightly varies with or without year_centered considered as RE, but story of signifigance remains the same

## testing the quality of the model
Jane_Site_Residuals2 <- simulateResiduals(Jane_Site_Model2, n = 1000)
plot(Jane_Site_Residuals2) #overall looks ok - good (but not perfect, similar to RE year_centered)
testDispersion(Jane_Site_Residuals2)
testZeroInflation(Jane_Site_Residuals2)

#lets use AIC testing here to compare all models and see if using the RE year_centered is viable
AIC(Jane_Site_Model, Jane_Site_Model2, Main_JaneRandom_Model)
# df      AIC
# Jane_Site_Model        9 103.6000: with year_centered RE is slightly better fit
# Jane_Site_Model2       7 104.4374
# Main_JaneRandom_Model 11 100.0350: acts as control but not much better fit than continuous year (good sign)

#confirm model fit
performance::check_singularity(Jane_Site_Model) #good

performance::check_model(Jane_Site_Model) #looks good, nothing too concerning

####model year as numeric & continuous based off main outplants####
# Create a continuous prediction sequence for every site (for modelling rates and to obtain summary stats (CI's))
Jane_Site_Pred <- JaneRandomYearCentered %>%
  dplyr::select(site, outplant) %>%
  distinct() %>%
  tidyr::crossing(year_num = seq(1988, 2022, length.out = 100)) %>%
  mutate(year_centered = (year_num - 1988) / 10) #previously (year_centered = (year_num - 2005) / 10) which made the error bars look a lot cleaner but it doesn't match our modeled time frame
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
  geom_point(data = JaneRandomYearCentered, aes(x = year_num, y = site_mean_density, colour = outplant), size = 2.5, alpha = 0.8) +
  facet_wrap(~site, ncol = 5, 
             labeller = labeller(site = \(x) str_to_title(str_replace_all(x, "_", " ")))) + # edit names
  scale_colour_manual(values = c("yes" = "#CD64B5", "no" = "grey60"),
                      labels = c("yes" = "Outplant", "no" = "Non-outplant")) +
  scale_fill_manual(values = c("yes" = "#CD64B5", "no" = "grey60"),
                    guide = "none") +
  scale_x_continuous(breaks = c(1988, 2007, 2022)) +
  scale_y_continuous(limits = c(-0.2, NA), expand = expansion(mult = c(0, 0.05))) + #set a -0.2 y so points dont get cut off
  labs(x = "Year", y = expression("Abalone density (" * m^-2 * ")"), colour = "") +
  theme_classic() +
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 9),
        legend.position = "right") +
  guides(color = guide_legend(reverse = FALSE)) #reverse legend to get yes on top
Jane_Site_Figure

#save fig
# #ggsave(
#   filename = "figures/MSc/Abalone_Jane_Site_Figure.png",
#   plot = Jane_Site_Figure,
#   width = 20,
#   height = 16,
#   units = "cm",
#   dpi = 900,
#   bg = "white")

####rate calculation####
#attempt slope figure 
#Explanation: Because we are modelling on the Log scale we have altered our relationship to match linear predictions from our exponential growth ones! All the log transformation does is change the absolute value between our relationship, not the proportion of change (this is why we can back transform (exponent) our estimates to get "real life" values of effect). So, all this is to say because we are modelling in "log land" we can get a fixed slope in log space using emtrends from our real model.

# Extract fixed effects
fixed_effects <- fixef(Jane_Site_Model)$cond

# Extract site-specific random effects
site_effects <- ranef(Jane_Site_Model)$cond$site

# Extract site-specific slopes
Site_Slopes <- tibble(
  site = rownames(site_effects),
  random_slope = site_effects[, "year_centered"]
) %>%
  left_join(
    JaneRandomYearCentered %>% distinct(site, outplant),
    by = "site"
  ) %>%
  mutate(
    slope = fixed_effects["year_centered"] +
      random_slope +
      if_else(
        outplant == "yes",
        fixed_effects["year_centered:outplantyes"],
        0
      )
  )

Site_Slopes <- Site_Slopes %>%
  mutate(
    percent_change_decade = 100 * (exp(slope) - 1),
    trend = case_when(
      slope > 0 ~ "Increasing",
      slope < 0 ~ "Decreasing",
      TRUE ~ "No change"
    )
  ) %>%
  arrange(slope)

ggplot(Site_Slopes,
       aes(x = percent_change_decade,
           y = reorder(site, percent_change_decade),
           colour = outplant)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_point(size = 3) +
  scale_colour_manual(
    values = c(
      "no" = "grey65",
      "yes" = "#CD64B5FF"
    ),
    labels = c("No", "Yes")
  ) +
  labs(
    x = "Estimated change in density per decade (%)",
    y = "Site",
    colour = "Received outplants?"
  ) +
  theme_classic()

#lets add CI's -> we must bootstrap to do so

get_site_slopes <- function(model, data) {
  
  # Extract fixed effects
  fixed_effects <- fixef(model)$cond
  
  # Extract site-specific random slopes
  site_effects <- ranef(model)$cond$site
  
  # Identify the non-reference outplant category
  non_reference <- levels(data$outplant)[2]
  
  # Identify the interaction coefficient automatically
  interaction_name <- names(fixed_effects)[
    grepl("year_centered", names(fixed_effects)) &
      grepl("outplant", names(fixed_effects)) &
      grepl(":", names(fixed_effects))
  ]
  
  # Extract interaction estimate
  interaction_effect <- unname(fixed_effects[interaction_name])
  
  # Calculate site-specific slopes
  tibble(
    site = rownames(site_effects),
    random_slope = site_effects[, "year_centered"]
  ) %>%
    left_join(
      data %>% distinct(site, outplant),
      by = "site"
    ) %>%
    mutate(
      slope = unname(fixed_effects["year_centered"]) +
        random_slope +
        if_else(
          outplant == non_reference,
          interaction_effect,
          0
        )
    )
}

#extract slopes from original model
Site_Slopes <- get_site_slopes(
  Jane_Site_Model,
  JaneRandomYearCentered)

Site_Slopes #why is the grappler random slope -ve? 

#set seed 
set.seed(123)

# Number of bootstrap simulations
B <- 1000

# Simulate new response datasets from the fitted model
Simulated_Data <- simulate(
  Jane_Site_Model,
  nsim = B)

Bootstrap_Slopes <- map_dfr(seq_len(B), function(i) {
  
  # Replace the original response with the simulated response
  Boot_Data <- JaneRandomYearCentered
  Boot_Data$site_mean_density <- Simulated_Data[[i]]
  
  # Refit the original model
  Boot_Model <- tryCatch(
    update(Jane_Site_Model, data = Boot_Data),
    error = function(e) NULL)
  
  # Skip failed model fits
  if (is.null(Boot_Model)) return(NULL)
  
  # Skip fits with a non-positive-definite Hessian
  if (!isTRUE(Boot_Model$sdr$pdHess)) return(NULL)
  
  # Extract the site-specific slopes
  get_site_slopes(Boot_Model, Boot_Data) %>%
    mutate(iteration = i)
  
})

# Calculate 95% CIs and convert slopes to annual percentage change
Site_Slope_CI <- Bootstrap_Slopes %>%
  group_by(site) %>%
  summarise(
    lower_slope = quantile(slope, 0.025, na.rm = TRUE),
    upper_slope = quantile(slope, 0.975, na.rm = TRUE),
    successful_bootstraps = n(),
    .groups = "drop"
  ) %>%
  left_join(Site_Slopes, by = "site") %>%
  mutate(
    # Convert decade slopes to annual slopes
    annual_slope = slope / 10,
    lower_annual = lower_slope / 10,
    upper_annual = upper_slope / 10,
    
    # Convert to percentage change per year
    percent_change = 100 * (exp(annual_slope) - 1),
    lower_percent = 100 * (exp(lower_annual) - 1),
    upper_percent = 100 * (exp(upper_annual) - 1)
  ) %>%
  arrange(annual_slope)

Site_Slope_CI

#plot
ggplot(
  Site_Slope_CI,
  aes(x = percent_change, y = reorder(site, percent_change),colour = outplant)) +
  geom_vline(xintercept = 0, linetype = "dashed",colour = "grey50") +
  geom_segment(aes(x = lower_percent, xend = upper_percent, y = reorder(site, percent_change), yend = reorder(site, percent_change)), linewidth = 0.8) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("no" = "grey65", "yes" = "#CD64B5FF"),
    labels = c("no" = "Non-outplant", "yes" = "Outplant")) +
  guides(color = guide_legend(reverse = FALSE)) + #reverse legend to get yes on top
  labs(x = "Estimated change in density per year (%)", y = "Site", colour = "") +
  scale_y_discrete(labels = \(x) str_to_title(str_replace_all(x, "_", " ")))+
  theme_classic()
#CI's are too large to detect which is unsurprising given the limitted temporal data, notice even if we shift to the medium outplant dataset we dont increase observations so this wont fix that...

#check boot straps of all sites (notice outplants have the only real discrepancy)
Bootstrap_Slopes <- Bootstrap_Slopes %>%
  group_by(site) %>%
  mutate(site_mean_slope = mean(slope, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(site = fct_reorder(site, site_mean_slope, .desc = TRUE)) #reorder names so that most +ve values appear first

Bootstrap_Slopes %>%
  ggplot(aes(x = slope, fill = outplant)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~site, scales = "free_y",
             labeller = labeller(site = \(x) str_to_title(str_replace_all(x, "_", " ")))) + # edit names) +
  theme_classic() +
  scale_fill_manual(values = c("yes" = "#CD64B5FF",  # Outplant sites
                               "no" = "grey65"),       # Non-outplant sites
                    labels = c("yes" = "Outplant", "no" = "Non-outplant"))+
  guides(fill = guide_legend(reverse = FALSE)) + #reverse legend to get yes first
  labs(x = "Site-specific temporal slope (log-density scale)", y = "Bootstrap density", fill = "")+
  theme(legend.position = "bottom")
#This fig explains why the CIs all overlap with 0. Sample size too small for the simulated data to generate smaller CIs

#now make a grouped forest plot rather than a by site comparison
Group_Slopes <- emtrends(
  Jane_Site_Model,
  specs = ~ outplant,
  var = "year_centered")

summary(Group_Slopes, infer = c(TRUE, TRUE))
Group_Slope_CI <- as.data.frame(confint(Group_Slopes)) #make df to call

#convert slope into change per year rather than decade
Group_Slope_Data <- Group_Slope_CI %>%
  mutate(
    annual_slope = year_centered.trend / 10,
    lower_annual = asymp.LCL / 10,
    upper_annual = asymp.UCL / 10,
    
    percent_change = 100 * (exp(annual_slope) - 1),
    lower_percent = 100 * (exp(lower_annual) - 1),
    upper_percent = 100 * (exp(upper_annual) - 1),
    
    outplant = factor(
      outplant,
      levels = c("no", "yes"),
      labels = c("Non-outplant", "Outplant")
    )
  )

Group_Slope_Data

#grouped forest plot
ggplot(
  Group_Slope_Data,
  aes(
    x = percent_change,
    y = outplant,
    colour = outplant)) + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = lower_percent, xend = upper_percent, yend = outplant), linewidth = 1) +
  geom_point(size = 4) +
  scale_colour_manual(values = c("Non-outplant" = "grey65", "Outplant" = "#CD64B5FF")) +
  labs(x = "Estimated rate of change in density per year (%)", y = NULL, colour = "Site group") +
  theme_classic() +
  theme(legend.position = "none")
