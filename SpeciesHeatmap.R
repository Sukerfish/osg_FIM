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
MonthlyMaster <- sqlFetch(conn, "hsdb_tbl_corp_physical_master") %>%
  filter_at(#remove non-FIM project data
    vars(starts_with("Project")), 
    any_vars(str_detect(., regex("AM",ignore_case = TRUE))))
odbcClose(conn)

# select gear and bay
RefsList <- MonthlyMaster %>%
  filter(Gear == 20) %>%
  filter(str_detect(Reference,"^TBM"))

# Establish Zone filter
ZoneFilter = c("A"
               ,"B"
               ,"C"
               ,"D"
               ,"E"
               )

# clean up biology and merge with selected references --------
CleanBio   <- osg_CleanBio(BioNumFull, RefsList)

# merge difficult to ID taxa and convert NODC code to human readable
CleanHRBio <- osg_ComBio(CleanBio, SpeciesList)


##### Tidy up for Analyses #####

#Red Tide Switch
#Coded as months of interest: e.g., Jul-Oct
RTS <- c(7:10)
StartYear <- 1998
EndYear <- 2017

# select only RT months and other filters
RT_Abund <- CleanHRBio %>%
  mutate(SAVsq   = BottomVegCover^2) %>%
  mutate(RTLogic = "A") %>%
  mutate(RTLogic = replace(RTLogic, year == 2005, "B")) %>%
  mutate(RTLogic = replace(RTLogic, year > 2005, "C")) %>%
  filter(month %in% RTS) %>%
  filter(#remove early years
    year >= StartYear) %>%
  filter(Zone %in% ZoneFilter)

# spread via scientific name
HaulFull <- RT_Abund %>%
  spread(Scientificname, N2) %>%
  replace(is.na(.),0)
# set zone as factor for rank abundance stuff
HaulFull$Zone      <- as.factor(HaulFull$Zone)
HaulFull$RTLogic   <- as.factor(HaulFull$RTLogic)
#HaulFull$RTLogic <- ordered(HaulFull$RTLogic, levels = c("Before", "During", "After"))
HaulFull$Stratum   <- as.factor(HaulFull$Stratum)
HaulFull$Reference <- as.character(HaulFull$Reference)

# ZONE FILTER #
#HaulFull <- HaulFull %>%
  #filter(Zone %in% c("D","E"))
  #filter(Stratum == "B")

# pull all abundance data out and
# REMOVE hauls with NO abundance data
HaulAbun <- HaulFull %>%
  subset(select = -c(month:RTLogic)) %>%
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric)))) %>%
  filter(total != 0) %>%
  subset(select = -c(total)) %>%
  as.data.frame()

HaulZeroAbun <- HaulFull %>%
  subset(select = -c(month:RTLogic)) %>%
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric)))) %>%
  filter(total == 0) %>%
  subset(select = -c(total)) %>%
  as.data.frame()

HaulFullClean <- HaulFull %>%
  inner_join(HaulAbun)

# pull all environmental data out
HaulEnv <- HaulFullClean %>%
  select(Reference:RTLogic) %>%
  as.data.frame()

library(vegan)

test <- RT_Abund %>%
  anti_join(HaulZeroAbun) %>%
  group_by(Reference) %>%
  summarise(Scientificname, Zone, N2)

  filter(Reference != c(HaulZeroAbun$Reference))
  anti_join(HaulZeroAbun)
  subset(select = -c(Reference))

ggplot(test, aes(Zone, Scientificname)) +
  geom_tile(aes(fill = N2), colour = "grey50")

