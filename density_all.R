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

# clean up biology and merge with selected references
CleanBio   <- osg_CleanBio(BioNumFull, RefsList)

# merge difficult to ID taxa and convert NODC code to human readable
CleanHRBio <- osg_ComBio(CleanBio, SpeciesList)

##### Tidy up for Analyses #####

#Red Tide Switch
#Coded as months of interest: e.g., Jul-Oct
RTS <- c(6:9)
#RTS <- c(1:3)
StartYear <- 1998
EndYear <- 2017
effort <- (140/100)

# select only RT months and other filters
RT_Abund <- CleanHRBio %>%
  mutate(season = "Winter") %>%
  mutate(season = replace(season, month %in% c(4:7), "Spring")) %>%
  #mutate(CPUE = N2/effort) %>%
  mutate(RTLogic = "Before") %>%
  mutate(RTLogic = replace(RTLogic, year == 2005, "During")) %>%
  mutate(RTLogic = replace(RTLogic, year > 2005, "After")) %>%
  #filter(month %in% RTS) %>%
  filter(#remove early years
    year >= StartYear) %>%
  filter(Zone %in% ZoneFilter)



# spread via scientific name
HaulFullTBM <- RT_Abund %>%
  spread(Scientificname, N2) %>%
  replace(is.na(.),0)
# set zone as factor for rank abundance stuff
HaulFullTBM$Zone      <- as.factor(HaulFullTBM$Zone)
HaulFullTBM$RTLogic   <- as.factor(HaulFullTBM$RTLogic)
#HaulFull$RTLogic <- ordered(HaulFull$RTLogic, levels = c("Before", "During", "After"))
HaulFullTBM$Stratum   <- as.factor(HaulFullTBM$Stratum)
HaulFullTBM$Reference <- as.character(HaulFullTBM$Reference)

HaulAbunTBM <- HaulFullTBM %>%
  subset(select = -c(month:RTLogic)) %>%
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric)))) %>%
  filter(total != 0) %>%
  #subset(select = -c(total)) %>%
  as.data.frame()

HaulZeroAbunTBM <- HaulFullTBM %>%
  subset(select = -c(month:RTLogic)) %>%
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric)))) %>%
  filter(total == 0) %>%
  #subset(select = -c(total)) %>%
  as.data.frame()

HaulFullCleanTBM <- HaulFullTBM %>%
  inner_join(HaulAbunTBM) %>%
  mutate(CPUE = total/effort)


# 
# # spread via scientific name and CPUE
# HaulFull <- RT_Abund %>%
#   spread(Scientificname, CPUE) %>%
#   replace(is.na(.),0)
# # set zone as factor for rank abundance stuff
# HaulFull$Zone    <- as.factor(HaulFull$Zone)
# HaulFull$RTLogic <- as.factor(HaulFull$RTLogic)
# HaulFull$RTLogic <- ordered(HaulFull$RTLogic, levels = c("Before", "During", "After"))
# HaulFull$Stratum <- as.factor(HaulFull$Stratum)
# 
# # # pull all environmental data out
# HaulEnv <- HaulFull %>%
#   select(Reference:RTLogic) %>%
#   as.data.frame()
# 
# # pull all abundance data out
# HaulAbun <- HaulFull %>%
#   subset(select = -c(Reference:RTLogic)) %>%
#   as.data.frame()

# summary information about the samples
HaulCount <- HaulFullCleanTBM %>%
  count(year)
HaulMin <- min(HaulCount$n)
#ZoneCount <- length(ZoneFilter)

library(vegan)
library(vegan3d)
library(goeveg)
library(scales)
library(ggrepel)

##### Yearly Resampled Species Richness ####

# Resample using the minimum samples in the paired values
HaulSub <- HaulFullCleanTBM %>%
  group_by(year) %>%
  # resample point
  sample_n(HaulMin, replace = FALSE) %>%
  ungroup()

# calculate average richness per haul using the subsampled data
# initialize
Haul_AvgRich <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("year", "avg_richness"))
for (i in 1:((EndYear-StartYear)+1)){
  temp_df <- HaulSub %>%
    # filter by year and remove extra columns
    filter(year == i+(StartYear-1)) %>%
    subset(select = CPUE)
    # assign year to row
    Haul_AvgRich[i,1]  <- i+(StartYear-1)
    # calculate richness per haul and average it before placement in matrix
    Haul_AvgRich[i,2] <- mean(temp_df$CPUE)
  }

# summary stats and plotting processing
Haul_tsRich <- Haul_AvgRich %>%
  # calculate as anomalies
  mutate(anom.Rich = avg_richness - mean(avg_richness, na.rm = TRUE))

RichnessAnom_Metrics <- Haul_tsRich %>%
  # calculate mean/CIs
  summarise(mean.Rich = mean(avg_richness, na.rm = TRUE),
            sd.anom.Rich = sd(anom.Rich, na.rm = TRUE),
            n.Rich = n()) %>%
  mutate(se.anom.Rich = sd.anom.Rich / sqrt(n.Rich),
         lower.ci.anom.Rich = 0 - (1.96 * se.anom.Rich),
         upper.ci.anom.Rich = 0 + (1.96 * se.anom.Rich))

####### Plots ######
# Richness over time plot
ggplot(Haul_tsRich, aes(x=year))+
  geom_line(aes(y=anom.Rich))+
  geom_hline(aes(yintercept = 0, color = "darkred"))+
  geom_ribbon(data=merge(RichnessAnom_Metrics, Haul_tsRich),
              aes(ymin=lower.ci.anom.Rich, 
                  ymax=upper.ci.anom.Rich), 
              linetype=2, alpha=0.1, color="purple")+
  theme_bw()+
  labs(title = "Animals per 100m^2 Over Time")+
  theme(axis.title = element_text(face = "bold", size = "12"), 
        axis.text = element_text(size = "12"),
        strip.text = element_text(face = "bold", size = "14"))+
  scale_x_continuous(name = "Year", 
                     limits = c(StartYear,EndYear),
                     breaks = seq(StartYear,EndYear, 2))+
  scale_y_continuous(name = "Average Density Anomaly from Long Term Mean",
                     limits = c(-200,250))+
  theme(legend.position="none")
  #facet_grid(Zone ~ .)