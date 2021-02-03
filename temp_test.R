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

# Establish Zone filters
APMZoneFilter = data.frame(system = "AP", 
                           Zone = c("A"
                                    ,"B"
                           ))

CKMZoneFilter = data.frame(system = "CK", 
                           Zone = c("B"
                                    ,"C"
                           ))

TBMZoneFilter = data.frame(system = "TB", 
                           Zone = c("A"
                                    ,"B"
                                    ,"C"
                                    ,"D"
                                    ,"E"
                           ))

CHMZoneFilter = data.frame(system = "CH", 
                           Zone = c("A"
                                    ,"B"
                                    ,"C"
                                    ,"D"
                           ))

ZoneFilter <- bind_rows(APMZoneFilter, 
                        CKMZoneFilter, 
                        TBMZoneFilter, 
                        CHMZoneFilter)


# clean up biology and merge with selected references
CleanBio   <- osg_CleanBio(BioNumFull, RefsList)

# merge difficult to ID taxa and convert NODC code to human readable
CleanHRBio <- osg_ComBio(CleanBio, SpeciesList)

FullHydro <- inner_join(HydroLab, RefsList, "Reference")

TidyHydro <- FullHydro %>%
  group_by(Reference) %>%
  summarise(depth = mean(Depth),
            temp  = mean(Temperature),
            cond  = mean(Conductivity),
            pH    = mean(pH),
            sal   = mean(Salinity),
            DO    = mean(DissolvedO2))

##### Tidy up for Analyses #####

#Red Tide Switch
#Coded as months of interest: e.g., Jul-Oct
#RTS <- c(6:9)
RTS <- c(1:3)
StartYear <- 2001
EndYear <- 2017
effort <- (140/100)

# select only RT months and other filters
RT_Abund <- CleanHRBio %>%
  mutate(system = if_else(str_detect(Reference, "^APM"), "AP",
                          if_else(str_detect(Reference, "^CKM"), "CK",
                                  if_else(str_detect(Reference, "^TBM"), "TB",
                                          "CH")))) %>%
  mutate(season = "Winter") %>%
  mutate(season = replace(season, month %in% c(4:7), "Spring")) %>%
  mutate(RTLogic = "Before") %>%
  mutate(RTLogic = replace(RTLogic, year == 2005, "During")) %>%
  mutate(RTLogic = replace(RTLogic, year > 2005, "After")) %>%
  filter(month %in% RTS) %>%
  filter(#remove early years
    year >= StartYear) %>%
  inner_join(ZoneFilter)

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

HaulAbun <- HaulFull %>%
  subset(select = -c(month:RTLogic)) %>%
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric)))) %>%
  filter(total != 0) %>%
  #subset(select = -c(total)) %>%
  as.data.frame()

HaulZeroAbun <- HaulFull %>%
  subset(select = -c(month:RTLogic)) %>%
  rowwise() %>%
  mutate(total = sum(c_across(where(is.numeric)))) %>%
  filter(total == 0) %>%
  #subset(select = -c(total)) %>%
  as.data.frame()

HaulFullClean <- HaulFull %>%
  inner_join(HaulAbun) %>%
  mutate(CPUE = total/effort) %>%
  inner_join(TidyHydro)

# summary information about the samples
HaulCount <- HaulFullClean %>%
  group_by(system) %>%
  count(year)
HaulMin <- HaulCount %>%
  summarise(minHaul = min(n))

#ZoneCount <- length(ZoneFilter)

library(vegan)
library(vegan3d)
library(goeveg)
library(scales)
library(ggrepel)

##### Yearly Resampled Species density ####

# Resample using the minimum samples in the paired values
HaulSub <- HaulFullClean %>%
  left_join(HaulMin) %>%
  group_by(system, year) %>%
  # resample point
  sample_n(minHaul, replace = FALSE) %>%
  subset(select = -c(minHaul)) %>%
  ungroup()

TempAnom <- HaulSub %>%
  subset(select = c(system, year, temp)) %>%
  group_by(system, year) %>%
  summarise(avg_temp = mean(temp, na.rm = TRUE),
            n_temp = n()) %>%
  ungroup() %>%
  group_by(system) %>%
  mutate(LTM_temp = mean(avg_temp),
         sd_temp = sd(avg_temp, na.rm = TRUE)) %>%
  mutate(anom_temp = avg_temp - LTM_temp,
         se_temp = sd_temp/sqrt(n_temp),
         lower.ci.anom.temp = 0 - (1.96 * se_temp),
         upper.ci.anom.temp = 0 + (1.96 * se_temp))

# # calculate average richness per haul using the subsampled data
# # initialize
# Haul_AvgRich <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), 
#                          c("system", "year", "avg_richness"))
# for (i in 1:length(syslist)){
#   temp_df <- HaulSub %>%
#     filter(system == syslist[,i])
# }
# for (i in 1:((EndYear-StartYear)+1)){
#   temp_df <- HaulSub %>%
#     # filter by year and remove extra columns
#     filter(year == i+(StartYear-1)) %>%
#     subset(select = temp)
#     # assign year to row
#     Haul_AvgRich[i,1]  <- i+(StartYear-1)
#     # calculate richness per haul and average it before placement in matrix
#     Haul_AvgRich[i,2] <- mean(temp_df$temp, na.rm = TRUE)
#   }
# 
# # summary stats and plotting processing
# Haul_tsRich <- Haul_AvgRich %>%
#   # calculate as anomalies
#   mutate(anom.Rich = avg_richness - mean(avg_richness, na.rm = TRUE))
# 
# RichnessAnom_Metrics <- Haul_tsRich %>%
#   # calculate mean/CIs
#   summarise(mean.Rich = mean(avg_richness, na.rm = TRUE),
#             sd.anom.Rich = sd(anom.Rich, na.rm = TRUE),
#             n.Rich = n()) %>%
#   mutate(se.anom.Rich = sd.anom.Rich / sqrt(n.Rich),
#          lower.ci.anom.Rich = 0 - (1.96 * se.anom.Rich),
#          upper.ci.anom.Rich = 0 + (1.96 * se.anom.Rich))

####### Plots ######
# Richness over time plot

ggplot(data=TempAnom,
       aes(x=year, y=anom_temp, color=system)) +
  geom_ribbon(
    aes(ymin=lower.ci.anom.temp, 
        ymax=upper.ci.anom.temp), 
    linetype=2, alpha=0.1, color="purple")+
  theme(legend.position="none") +
  facet_wrap(as.factor(TempAnom$system), scales = "free") +
  geom_line()
# 
# ggplot(testtt, aes(x=year))+
#   geom_line(aes(y=anom_temp))+
#   geom_hline(aes(yintercept = 0, color = "darkred"))+
#   geom_ribbon(
#               aes(ymin=lower.ci.anom.temp, 
#                   ymax=upper.ci.anom.temp), 
#               linetype=2, alpha=0.1, color="purple")+
#   theme_bw()+
#   labs(title = "Average haul water temperature over time")+
#   theme(axis.title = element_text(face = "bold", size = "12"), 
#         axis.text = element_text(size = "12"),
#         strip.text = element_text(face = "bold", size = "14"))+
#   scale_x_continuous(limits = c(StartYear,EndYear),
#                      breaks = seq(StartYear,EndYear, 2))+
#   scale_y_continuous(name = "Temperature Anomaly from Long Term Mean",
#                      limits = c(-4,4))+
#   theme(legend.position="none") +
#   facet_grid(system)
# 
