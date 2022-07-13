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
HaulFull <- HaulFull %>%
  filter(Zone %in% c("D","E"))
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

HaulFullClean <- HaulFull %>%
  inner_join(HaulAbun)

# pull all environmental data out
HaulEnv <- HaulFullClean %>%
  select(Reference:RTLogic) %>%
  as.data.frame()

# # summary information about the samples using HaulAbun
# HaulCount <- HaulFullClean %>%
#   count(year, Zone)
# HaulMin <- min(HaulCount$n)
# ZoneCount <- length(ZoneFilter)

# # reorganize to include richness variable
# HaulRich <- HaulFull %>%
#   subset(select = -c(month:RTLogic))
# HaulRich$Richness <- specnumber(HaulRich, HaulRich$Reference)
# HaulRich <- HaulRich %>%
#   subset(select = c(Reference, Richness)) %>%
#   merge(HaulEnv)
# 
# # reorganize to count all organisms captured per haul
# HaulTotAbun <- RT_Abund %>%
#   group_by(Reference) %>%
#   summarise(totalAb = sum(N2)) %>%
#   merge(HaulEnv)

library(vegan)
#library(vegan3d)
#library(goeveg)
#library(scales)
#library(ggrepel)
#library(MASS)
#library(nlme)

##### Mixed Effects models ####

# # Resample using the minimum samples in the paired values
# HaulSub <- HaulFull %>%
#   group_by(year, Zone) %>%
#   # resample point
#   sample_n(HaulMin, replace = FALSE) %>%
#   ungroup()


test <- HaulAbun %>%
  subset(select = -c(Reference))

#test <- decostand(test,"pa")

example_NMDS=metaMDS(test, distance = "bray", k = 2)
adon.results = adonis(test ~ HaulEnv$RTLogic, method="bray",perm=999)
print(adon.results)

stressplot(example_NMDS)
plot(example_NMDS)
ordiplot(example_NMDS,type="n")
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",cex=1.25,air=0.01)

treat=c(rep("Treatment1",5),rep("Treatment2",5))
ordiplot(example_NMDS,type="n")
ordihull(example_NMDS,groups=treat,draw="polygon",col="grey90",label=F)
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",col=c(rep("green",5),rep("blue",5)),
         air=0.01,cex=1.25)


# Resample using the minimum samples in the paired values
HaulSub <- HaulFullClean %>%
  #inner_join(HaulRich) %>%
  #inner_join(HaulTotAbun) %>%
  group_by(year, Zone) %>%
  # resample point
  sample_n(HaulMin, replace = FALSE) %>%
  ungroup()

HaulSubSppAvgFull <- HaulSub %>%
  group_by(year, Zone) %>%
  summarise(
    across(where(is.numeric), mean)) %>%
  select(-c(month:SAVsq)) %>%
  #select(-c(Richness:totalAb)) %>%
  ungroup()
  
HaulSubEnvAvg <- HaulFull %>%
  select(c(Reference:RTLogic))

HaulSubSppAvg <- HaulFull %>%
  select(-c(Reference:RTLogic))

HaulSubSppAvg[HaulSubSppAvg >0] <- 1