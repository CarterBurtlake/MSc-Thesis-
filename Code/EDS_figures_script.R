#Figures for EDS report

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
  scale_y_continuous(limits = c(-1, 3), breaks = c(0,0.5,1,1.5,2,2.5, 3))+ #play with this so that we can see the full range of CI but also have real and relevant y axis limits (without we plot into -ve)
  facet_wrap(~site) #show by site

#now we need to pull the slope value to get an average change in abalone density group

#Fig 1b) --------------------------------------------------------------------------------------------

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
  ungroup()

#Now lets plot density through time per site for the combined data set
ggplot(data = Fig1_JW_RLS, aes(x = as.numeric(year), y = site_mean_density, colour = outplant)) +
  geom_point()+
  geom_smooth(method = "lm") +
  labs(x = "Year", y = "Abalone density (" ~m^-2*")", colour = "Recieved outplants?") + #add titles
  theme_classic()+ #white background
  guides(color = guide_legend(reverse = TRUE))+ #put yes status on the top of the legend
  scale_y_continuous(limits = c(-1, 3), breaks = c(0,0.5,1,1.5,2,2.5, 3))+ #play with this so that we can see the full range of CI but also have real and relevant y axis limits (without we plot into -ve)
  facet_wrap(~site) #show by site
