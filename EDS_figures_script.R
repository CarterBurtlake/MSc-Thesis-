#Figures for EDS report

#Figure 1 -------------------------------------------------------------------------------------------

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
  facet_wrap(~site)

