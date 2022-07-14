#### Import and Filter Data #####
 
library(tidyverse)
library(lubridate)

load('GearCode20Refresh.Rdata')

#establish years of interest
YearFilter <- ZoneFilter %>%
  select(system, StartYear, EndYear) %>%
  unique()

#establish seasons of interest
SOI <- c(
  "summer",
  "winter"
)

#expand sampling date and attach years of interest
HydroList <- TidyHydro %>%
  mutate(month = month(Sampling_Date), year = year(Sampling_Date)) %>%
  left_join(YearFilter) %>%
  #label specific months as seasons
  mutate(season = if_else(month %in% c(4:5), "spring",
                          if_else(month %in% c(6:9), "summer",
                                  if_else(month %in% c(10:11), "fall",
                                          "winter")))) %>%
  #tag everything with year associated with sampling season, with December
  #as the preceding year's data
  mutate(seasonYear = if_else(month %in% c(12), year + 1, year)) %>%
  #remove the incomplete final seasonYear of winter data
  filter(season != "winter" | seasonYear != EndYear + 1) %>%
  #filter by years of interest (all begin before, so first year of interest is
  #not truncated)
  filter(seasonYear >= StartYear) %>%
  #remove year of interest logic
  subset(select = -c(StartYear, EndYear)) %>%
  #filter by seasons of interest
  filter(season %in% SOI) 

#Zone filtering setup

#establish zones of interest by estuary
ZoneLogic <- ZoneFilter %>%
  subset(select = c(system, Zone)) %>%
  mutate(systemZone = str_c(system, "_", Zone)) %>%
  subset(select = systemZone)

#build new reference/Zone dataframe for processing
TidyRefsList <- RefsList %>%
  subset(select = c(Reference, Zone)) %>%
  inner_join(HydroList) %>%
  #sanitize Zone entries
  mutate(Zone = str_trim(Zone, side  = "both")) %>%
  mutate(Zone = str_to_upper(Zone)) %>%
  mutate(systemZone = str_c(system, "_", Zone)) %>%
  #filter by system_zones in logic vector
  filter(systemZone %in% ZoneLogic$systemZone) %>%
  subset(select = c(Reference, system, season, seasonYear, systemZone)) 

TidyBio <- TidyRefsList %>%
  left_join(CleanHRBio, by = "Reference")

save(CleanHRBio, 
     RefsList, 
     HydroList, 
     ZoneLogic, 
     TidyRefsList, 
     TidyBio, 
     file = "TidyGearCode20.Rdata")
