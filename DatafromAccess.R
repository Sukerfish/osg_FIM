#### Import and Filter Data #####

library(tidyverse)
library(lubridate)
library(RODBC)
library(devtools)
source_url("https://github.com/Sukerfish/osg_FIM/blob/master/Functions/osg_CleanBio.R?raw=TRUE")
source_url("https://github.com/Sukerfish/osg_FIM/blob/master/Functions/osg_ComBio.R?raw=TRUE")

conn <- odbcConnect("FIM_db")
SpeciesList   <- sqlFetch(conn, "hsdb_tbl_corp_ref_species_list")
SpeciesList$NODCCODE <- as.character(SpeciesList$NODCCODE)
BioNumFull    <- sqlFetch(conn, "hsdb_tbl_corp_biology_number")
HydroLab      <- sqlFetch(conn, "hsdb_tbl_corp_hydrolab")
MonthlyMaster <- sqlFetch(conn, "hsdb_tbl_corp_physical_master") %>%
  filter_at(#remove non-FIM project data
    vars(starts_with("Project")), 
    any_vars(str_detect(., regex("AM",ignore_case = TRUE))))
odbcClose(conn)

# select gear and bays
RefsList <- MonthlyMaster %>%
  filter(Gear == 20) %>%
  filter(str_detect(Reference, paste(c(
    "^APM",
    "^CKM",
    "^TBM",
    "^CHM"
  ), collapse = '|')))

# Establish Zone and Year filters
APMZoneFilter = data.frame(system = "AP", 
                           Zone = c("A"
                                    ,"B"
                           ),
                           StartYear = 2001,
                           EndYear = 2017)

CKMZoneFilter = data.frame(system = "CK", 
                           Zone = c("B"
                                    ,"C"
                           ),
                           StartYear = 2001,
                           EndYear = 2017)

TBMZoneFilter = data.frame(system = "TB", 
                           Zone = c("A"
                                    ,"B"
                                    ,"C"
                                    ,"D"
                                    ,"E"
                           ),
                           StartYear = 1998,
                           EndYear = 2017)

CHMZoneFilter = data.frame(system = "CH", 
                           Zone = c("A"
                                    ,"B"
                                    ,"C"
                                    ,"D"
                           ),
                           StartYear = 1998,
                           EndYear = 2017)

ZoneFilter <- bind_rows(APMZoneFilter, 
                        CKMZoneFilter, 
                        TBMZoneFilter, 
                        CHMZoneFilter)


# clean up biology and merge with selected references
CleanBio   <- osg_CleanBio(BioNumFull, RefsList)

# merge difficult to ID taxa and convert NODC code to human readable
CleanHRBio <- osg_ComBio(CleanBio, SpeciesList)

# pull in Hydro data for all Refs of Interest
FullHydro <- inner_join(HydroLab, RefsList, "Reference")

#summarize Hydro data for each Reference
TidyHydro <- FullHydro %>%
  select(c('Reference',
  'Sampling_Date',
  'Depth',
  'Temperature',
  'Conductivity',
  'pH',
  'Salinity',
  'DissolvedO2')) %>%
  group_by(Reference, Sampling_Date) %>%
  summarise(across(, mean)) %>%
  mutate(system = if_else(str_detect(Reference, "^APM"), "AP",
                          if_else(str_detect(Reference, "^CKM"), "CK",
                                  if_else(str_detect(Reference, "^TBM"), "TB",
                                          "CH"))))

#save(CleanHRBio, RefsList, TidyHydro, ZoneFilter, file = "GearCode20.Rdata")