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

#create logic code for zone filtering
# ZoneLogic <- ZoneFilter %>%
#   subset(select = c(system, Zone)) %>%
#   group_by(system) %>%
#   summarise(code=paste(Zone,collapse=", "))

TidyRefsList <- RefsList %>%
  subset(select = c(Reference, Zone)) %>%
  inner_join(HydroList) %>%
  #inner_join(ZoneLogic) %>%
  mutate(Zone = str_trim(Zone, side  = "both")) %>%
  filter(Zone %in% ZoneFilter$Zone) %>%
  subset(select = c(Reference, system, season, seasonYear, Zone))

TidyBio <- TidyRefsList %>%
  left_join(CleanHRBio, by = "Reference")

save(CleanHRBio, 
     RefsList, 
     HydroList, 
     ZoneFilter, 
     TidyRefsList, 
     TidyBio, 
     file = "TidyGearCode20.Rdata")
