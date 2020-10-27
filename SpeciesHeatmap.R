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
TBMList <- MonthlyMaster %>%
  filter(Gear == 20) %>%
  filter(str_detect(Reference,"^TBM"))

# Establish Zone filter
ZoneFilterTBM = c("A"
               ,"B"
               ,"C"
               ,"D"
               ,"E"
               )

# select gear and bay
CHMList <- MonthlyMaster %>%
  filter(Gear == 20) %>%
  filter(str_detect(Reference,"^CHM"))

# Establish Zone filter
ZoneFilterCHM = c("A"
               ,"B"
               ,"C"
               ,"D"
               #,"E"
)


# clean up biology and merge with selected references --------
CleanBioTBM   <- osg_CleanBio(BioNumFull, TBMList)
CleanBioCHM   <- osg_CleanBio(BioNumFull, CHMList)

# merge difficult to ID taxa and convert NODC code to human readable
CleanHRBioTBM <- osg_ComBio(CleanBioTBM, SpeciesList)
CleanHRBioCHM <- osg_ComBio(CleanBioCHM, SpeciesList)


##### Tidy up for Analyses #####

#Red Tide Switch
#Coded as months of interest: e.g., Jul-Oct
RTS <- c(6:9)
StartYear <- 1998
EndYear <- 2017

# select only RT months and other filters
RT_AbundTBM <- CleanHRBioTBM %>%
  mutate(SAVsq   = BottomVegCover^2) %>%
  mutate(RTLogic = "A") %>%
  mutate(RTLogic = replace(RTLogic, year == 2005, "B")) %>%
  mutate(RTLogic = replace(RTLogic, year > 2005, "C")) %>%
  filter(month %in% RTS) %>%
  filter(#remove early years
    year >= StartYear) %>%
  filter(Zone %in% ZoneFilterTBM)

# select only RT months and other filters
RT_AbundCHM <- CleanHRBioCHM %>%
  mutate(SAVsq   = BottomVegCover^2) %>%
  mutate(RTLogic = "A") %>%
  mutate(RTLogic = replace(RTLogic, year == 2005, "B")) %>%
  mutate(RTLogic = replace(RTLogic, year > 2005, "C")) %>%
  filter(month %in% RTS) %>%
  filter(#remove early years
    year >= StartYear) %>%
  filter(Zone %in% ZoneFilterCHM)

# spread via scientific name
HaulFullTBM <- RT_AbundTBM %>%
  spread(Scientificname, N2) %>%
  replace(is.na(.),0)
# set zone as factor for rank abundance stuff
HaulFullTBM$Zone      <- as.factor(HaulFullTBM$Zone)
HaulFullTBM$RTLogic   <- as.factor(HaulFullTBM$RTLogic)
#HaulFull$RTLogic <- ordered(HaulFull$RTLogic, levels = c("Before", "During", "After"))
HaulFullTBM$Stratum   <- as.factor(HaulFullTBM$Stratum)
HaulFullTBM$Reference <- as.character(HaulFullTBM$Reference)

# spread via scientific name
HaulFullCHM <- RT_AbundCHM %>%
  spread(Scientificname, N2) %>%
  replace(is.na(.),0)
# set zone as factor for rank abundance stuff
HaulFullCHM$Zone      <- as.factor(HaulFullCHM$Zone)
HaulFullCHM$RTLogic   <- as.factor(HaulFullCHM$RTLogic)
#HaulFull$RTLogic <- ordered(HaulFull$RTLogic, levels = c("Before", "During", "After"))
HaulFullCHM$Stratum   <- as.factor(HaulFullCHM$Stratum)
HaulFullCHM$Reference <- as.character(HaulFullCHM$Reference)

# ZONE FILTER #
#HaulFull <- HaulFull %>%
  #filter(Zone %in% c("D","E"))
  #filter(Stratum == "B")

# pull all abundance data out and
# REMOVE hauls with NO abundance data
HaulAbunTBM <- HaulFullTBM %>%
  subset(select = -c(month:RTLogic)) %>%
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric)))) %>%
  filter(total != 0) %>%
  subset(select = -c(total)) %>%
  as.data.frame()

HaulZeroAbunTBM <- HaulFullTBM %>%
  subset(select = -c(month:RTLogic)) %>%
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric)))) %>%
  filter(total == 0) %>%
  subset(select = -c(total)) %>%
  as.data.frame()

HaulFullCleanTBM <- HaulFullTBM %>%
  inner_join(HaulAbunTBM)

# pull all abundance data out and
# REMOVE hauls with NO abundance data
HaulAbunCHM <- HaulFullCHM %>%
  subset(select = -c(month:RTLogic)) %>%
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric)))) %>%
  filter(total != 0) %>%
  subset(select = -c(total)) %>%
  as.data.frame()

HaulZeroAbunCHM <- HaulFullCHM %>%
  subset(select = -c(month:RTLogic)) %>%
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric)))) %>%
  filter(total == 0) %>%
  subset(select = -c(total)) %>%
  as.data.frame()

HaulFullCleanCHM <- HaulFullCHM %>%
  inner_join(HaulAbunCHM)

# pull all environmental data out
HaulEnvTBM <- HaulFullCleanTBM %>%
  select(Reference:RTLogic) %>%
  mutate(System = "TB") %>%
  as.data.frame()

# pull all environmental data out
HaulEnvCHM <- HaulFullCleanCHM %>%
  select(Reference:RTLogic) %>%
  mutate(System = "CH") %>%
  as.data.frame()

#library(vegan)

# testtwoTBM <- HaulAbunTBM %>%
#   subset(select=-c(Reference)) %>%
#   summarise_all(list(sum))
# 
# testtwoTBM <- sort(testtwoTBM, decreasing = TRUE)
# 
# testtwoCHM <- HaulAbunCHM %>%
#   subset(select=-c(Reference)) %>%
#   summarise_all(list(sum))
# 
# testtwoCHM <- sort(testtwoCHM, decreasing = TRUE)
#   
# test <- full_join(HaulAbunTBM, HaulAbunCHM) %>%
#   subset(select=-c(Reference)) %>%
#   summarise_all(list(sum))
# 
# # removes species not found in CHM currently
# test <- sort(test, decreasing = TRUE)
# 
# testenv <- full_join(HaulEnvCHM, HaulEnvTBM)

rankTBM <- RT_AbundTBM %>%
  anti_join(HaulZeroAbunTBM) %>%
  group_by(Scientificname) %>%
  summarise_at(vars(N2), sum) %>%
  mutate(rank = row_number(desc(N2))) %>%
  arrange(N2) %>%
  ungroup()
  
testTBM <- rankTBM %>%
  select(-N2) %>%
  inner_join(RT_AbundTBM) %>%
  group_by(rank, Zone, Scientificname) %>%
  summarise_at(vars(N2), sum) %>%
  mutate(pa = replace(N2, N2 > 0, 1)) %>%
  mutate(System = "TB") %>%
  mutate(Group = 3) %>%
  mutate(Group = replace(Group, rank < 101, 2)) %>%
  mutate(Group = replace(Group, rank < 51, 1))

testTBM$Scientificname <- factor(testTBM$Scientificname, levels = rankTBM$Scientificname)

#Same for CHM
rankCHM <- RT_AbundCHM %>%
  anti_join(HaulZeroAbunCHM) %>%
  group_by(Scientificname) %>%
  summarise_at(vars(N2), sum) %>%
  mutate(rank = row_number(desc(N2))) %>%
  arrange(N2) %>%
  ungroup()

testCHM <- rankCHM %>%
  select(-N2) %>%
  inner_join(RT_AbundCHM) %>%
  group_by(rank, Zone, Scientificname) %>%
  summarise_at(vars(N2), sum) %>%
  mutate(pa = replace(N2, N2 > 0, 1)) %>%
  mutate(System = "CH") %>%
  mutate(Group = 3) %>%
  mutate(Group = replace(Group, rank < 101, 2)) %>%
  mutate(Group = replace(Group, rank < 51, 1))

testCHM$Scientificname <- factor(testCHM$Scientificname, levels = rankCHM$Scientificname)

# testCHM <- RT_AbundCHM %>%
#   anti_join(HaulZeroAbunCHM) %>%
#   #filter(Scientificname %in% colnames(testtwoCHM[1:50])) %>%
#   #mutate(N2 = replace(N2, N2 > 0, 1)) %>%
#   group_by(Zone, Scientificname) %>%
#   summarise_at(vars(N2), sum) %>%
#   ungroup() %>%
#   group_by(Scientificname) %>%
#   #arrange(desc(N2), .by_group = TRUE) %>%
#   mutate(N2Sort = sum(N2)) %>%
#   ungroup() %>%
#   mutate(rank = min_rank(desc(N2Sort))) %>%
#   arrange(desc(N2Sort)) %>%
#   mutate(System = "CH")

# testBoth <- full_join(testTBM, testCHM)
# 
# ggplot(testBoth, aes(Zone, factor(Scientificname, levels = rank))) +
#   geom_tile(aes(fill = N2), colour = "grey50") +
#   guides(fill=FALSE) +
#   facet_grid(. ~ System)

ggplot(testTBM, aes(Zone, Scientificname)) +
  theme_bw()+
  theme(strip.text = element_blank(), 
        axis.text.x = element_text(vjust = 0.5, size = 12, face = "bold"), 
        axis.text.y = element_text(size = 11, face = "italic"), 
        axis.title = element_text(size = 13, face = "bold"),
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Tampa Bay Taxa Presence", x = "Zone",y = "Scientific Name")+
  geom_tile(aes(fill = pa), colour = "grey50") +
  guides(fill=FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~Group, scale = "free")

ggplot(testCHM, aes(Zone, Scientificname)) +
  theme_bw()+
  theme(strip.text = element_blank(), 
        axis.text.x = element_text(vjust = 0.5, size = 12, face = "bold"), 
        axis.text.y = element_text(size = 11, face = "italic"), 
        axis.title = element_text(size = 13, face = "bold"),
        plot.title = element_text(hjust = 0.5))+
  labs(title = "Charlotte Harbor Taxa Presence", x = "Zone",y = "Scientific Name")+
  geom_tile(aes(fill = pa), colour = "grey50") +
  guides(fill=FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~Group, scale = "free")

