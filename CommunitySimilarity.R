# Bray-Curtis similarity pairwise among all years with respect to first year of data
# Using average annual density of each taxa within each estuary and season
# 
# Model change in abundance through time for each taxa
# 
# Calculate coefficient of the slope of each
# 
# Potential coefficient distributions

library(tidyverse)

load('TidyGearCode20.Rdata')

CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, systemZone, Scientificname, N2)) %>%
  #filter(season == "summer" | season == "winter")
  filter(season == "winter")
  #filter(season == "summer")

SiteXSpeciesFull <- CleanHauls %>%
  group_by(Reference) %>%
  spread(Scientificname,N2) %>%
  ungroup() %>%
  subset(select = -c(Reference)) %>%
  group_by(system, season, seasonYear, systemZone) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

