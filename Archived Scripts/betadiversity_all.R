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
  group_by(Reference) %>%
  summarise(depth = mean(Depth),
            temp  = mean(Temperature),
            cond  = mean(Conductivity),
            pH    = mean(pH),
            sal   = mean(Salinity),
            DO    = mean(DissolvedO2))

##### Tidy up for Analyses #####

#Months of Interest
#Coded as months of interest: e.g., Jul-Oct

#MOI <- c(6:9)
MOI <- c(1:3)

effort <- (140/100)

# select only MOI and other filters
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
  filter(month %in% MOI) %>%
  #filter(#remove early years
    #year >= StartYear) %>%
  inner_join(ZoneFilter) %>%
  filter(year >= StartYear) %>%
  subset(select = -c(StartYear, EndYear))

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
#HaulFull$year      <- as.factor(HaulFull$year)

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

##### Yearly Resampling ####

# Resample using the minimum samples in the paired values
HaulSub <- HaulFullClean %>%
  left_join(HaulMin) %>%
  group_by(system, year) %>%
  # resample point
  sample_n(minHaul, replace = FALSE) %>%
  subset(select = -c(minHaul)) %>%
  ungroup()

HaulSub_Abun <- HaulSub %>%
  subset(select = -c(Stratum:ShoreDistance,
                     Reference:month,
                     season:RTLogic,
                     total:DO))

HaulSub_Env <- HaulSub %>%
  subset(select = c(Reference:RTLogic, 
                    total:DO))

##### Beta Diversity ########
library(codyn)
#library(betapart)

HaulSub_AvgAbun <- HaulSub_Abun %>%
  group_by(system, year) %>%
  summarise(across(everything(), mean)) %>%
  ungroup() %>%
  #select(which(!colSums(.[sapply(., is.numeric)], na.rm = TRUE) %in% 0)) %>%
  pivot_longer(!c(year, system), names_to = "species", values_to = "AvgAbun")
  
test <- HaulSub_AvgAbun %>%
  #filter(system == "TB") %>%
  select(-year) %>%
  group_by(system, species) %>% 
  summarise(across(everything(), mean)) %>%
  filter(AvgAbun != 0) %>%
  select(-AvgAbun) %>%
  inner_join(HaulSub_AvgAbun)
  #filter(system == "TB") %>%
  #select(-system)

# test <- HaulSub_AvgAbun %>%
#   filter(system == "AP") %>%
#   select(-system, -year) %>%
#   group_by(species) %>% 
#   summarise(across(everything(), mean)) %>%
#   filter(AvgAbun != 0) %>%
#   select(-AvgAbun) %>%
#   inner_join(HaulSub_AvgAbun) %>%
#   filter(system == "AP") %>%
#   select(-system)


total_turnover <- turnover(df = test, 
                         time.var = "year", 
                         species.var = "species", 
                         abundance.var = "AvgAbun",
                         replicate.var = "system"
                         )

rankshift <- rank_shift(df = test, 
                           time.var = "year", 
                           species.var = "species", 
                           abundance.var = "AvgAbun",
                           replicate.var = "system"
                      )

rankshift$year <- as.numeric(substr(rankshift$year_pair, 6, 9))

ggplot(rankshift, aes(year, MRS)) + 
  facet_wrap(~system) +
  geom_line(size= 1) + theme_bw() 

TurnAnom <- total_turnover %>%
  group_by(system) %>%
  mutate(LTM_total = mean(total),
            n_total = n(),
            sd_total = sd(total)) %>%
  mutate(se_total = sd_total/sqrt(n_total),
         lower.ci.anom.total = LTM_total - (se_total),
         upper.ci.anom.total = LTM_total + (se_total))

# RichAnom <- RichFull %>%
#   group_by(system, year) %>%
#   summarise(avg_rich = mean(richness),
#             n_rich = n()) %>%
#   ungroup() %>%
#   group_by(system) %>%
#   mutate(LTM_rich = mean(avg_rich),
#          sd_rich = sd(avg_rich)) %>%
#   mutate(anom_rich = avg_rich - LTM_rich,
#          se_rich = sd_rich/sqrt(n_rich),
#          lower.ci.anom.rich = 0 - (1.96 * se_rich),
#          upper.ci.anom.rich = 0 + (1.96 * se_rich))

app_turnover <- turnover(df = test, 
                           time.var = "year", 
                           species.var = "species", 
                           abundance.var = "AvgAbun",
                           replicate.var = "system",
                           metric = "appearance")

AppDis <- turnover(df = test, 
                   time.var = "year", 
                   species.var = "species", 
                   abundance.var = "AvgAbun",
                   replicate.var = "system",
                   metric = "disappearance") %>%
  left_join(app_turnover) %>%
  pivot_longer(!c(year, system), names_to = "direction", values_to = "prop")

AllRate <- rate_change_interval(df = test, 
                         time.var = "year", 
                         species.var = "species", 
                         abundance.var = "AvgAbun",
                         replicate.var = "system")
  
# 
# ###### Density Calculations ########
# HaulSub_Dens <- HaulSub_Abun %>%
#   mutate(across(!c(year, system), ~{.x/effort}))
# 
# DensFull <- HaulSub_Dens %>%
#   group_by(system, year) %>%
#   summarise(across(everything(), mean)) %>%
#   pivot_longer(!c(year, system), names_to = "species", values_to = "CPUE")
# 
# DensDelta <- ZoneFilter %>%
#   subset(select = -Zone) %>%
#   distinct() %>%
#   inner_join(DensFull) %>%
#   group_by(system) %>%
#   filter(year == StartYear |
#            year == StartYear +1 |
#            year == StartYear +2 |
#            year == EndYear |
#            year == EndYear -1 |
#            year == EndYear -2 ) %>%
#   mutate(time = ifelse(year == StartYear |
#                          year == StartYear +1 |
#                          year == StartYear +2 , "start", "end")) %>%
#   subset(select = -c(year, StartYear, EndYear)) %>%
#   ungroup() %>%
#   group_by(system, species, time) %>%
#   summarise(across(everything(), mean)) %>%
#   pivot_wider(names_from = c(time), 
#               values_from = c(CPUE)) %>%
#   mutate(delta = end - start) %>%
#   ungroup() %>%
#   mutate(chk = abs(end) + abs(start)) %>%
#   filter(chk != 0) %>%
#   subset(select = -c(chk)) %>%
#   arrange(desc(delta))
# 
# DensDeltaInf <- DensDelta %>%
#   filter(start == 0 | end == 0) %>%
#   group_by(system) %>%
#   arrange(desc(delta))
# 
# DensDeltaClean <- DensDelta %>%
#   filter(start > 0 & end > 0) %>%
#   mutate(deltaP = (end - start)/start) %>%
#   arrange(desc(deltaP))
#   
# ##### Temperature Anomaly ####
# TempFull <- HaulSub %>%
#   subset(select = c(system, year, temp))
# 
# TempAnom <- TempFull %>%
#   group_by(system, year) %>%
#   summarise(avg_temp = mean(temp, na.rm = TRUE),
#             n_temp = n()) %>%
#   ungroup() %>%
#   group_by(system) %>%
#   mutate(LTM_temp = mean(avg_temp),
#          sd_temp = sd(avg_temp, na.rm = TRUE)) %>%
#   mutate(anom_temp = avg_temp - LTM_temp,
#          se_temp = sd_temp/sqrt(n_temp),
#          lower.ci.anom.temp = 0 - (1.96 * se_temp),
#          upper.ci.anom.temp = 0 + (1.96 * se_temp))
# 
# ##### Richness Anomaly ######
# 
# RichFull <- HaulSub %>%
#   subset(select = -c(Stratum:ShoreDistance)) %>%
#   subset(select = -c(month)) %>%
#   subset(select = -c(season:RTLogic)) %>%
#   subset(select = -c(total:DO)) %>%
#   group_by(system, year) %>%
#   # mutate all non-zero things to be 1
#   # ignore Reference column, system and year used as grouping
#   mutate(across(-c(Reference), ~1 * (. >0))) %>%
#   ungroup() %>%
#   rowwise() %>%
#   mutate(richness = sum(c_across(where(is.numeric)))) %>%
#   ungroup() %>%
#   subset(select = c(Reference, system, year, richness))
# 
# RichDelta <- ZoneFilter %>%
#   subset(select = -Zone) %>%
#   distinct() %>%
#   inner_join(RichFull) %>%
#   subset(select = -c(Reference)) %>%
#   group_by(system) %>%
#   filter(year == StartYear |
#            year == StartYear +1 |
#            year == StartYear +2 |
#            year == EndYear |
#            year == EndYear -1 |
#            year == EndYear -2 ) %>%
#   mutate(time = ifelse(year == StartYear |
#                          year == StartYear +1 |
#                          year == StartYear +2 , "start", "end")) %>%
#   subset(select = -c(year, StartYear, EndYear)) %>%
#   ungroup() %>%
#   group_by(system, time) %>%
#   summarise(avg_rich = mean(richness))
# 
# RichAnom <- RichFull %>%
#   group_by(system, year) %>%
#   summarise(avg_rich = mean(richness),
#             n_rich = n()) %>%
#   ungroup() %>%
#   group_by(system) %>%
#   mutate(LTM_rich = mean(avg_rich),
#          sd_rich = sd(avg_rich)) %>%
#   mutate(anom_rich = avg_rich - LTM_rich,
#          se_rich = sd_rich/sqrt(n_rich),
#          lower.ci.anom.rich = 0 - (1.96 * se_rich),
#          upper.ci.anom.rich = 0 + (1.96 * se_rich))

####### Plots ######

# Temp anomaly over time
# ggplot(data=TempAnom,
#        aes(x=year, y=anom_temp, group=system, color=system)) +
#   geom_ribbon(
#     aes(ymin=lower.ci.anom.temp,
#         ymax=upper.ci.anom.temp),
#     linetype=2, alpha=0.1, color="purple")+
#   theme(legend.position="none") +
#   facet_wrap(as.factor(TempAnom$system), scales = "free") +
#   geom_line()
# 
# #Richness anomaly over time
# ggplot(data=RichAnom,
#        aes(x=year, y=anom_rich, group=system, color=system)) +
#   geom_ribbon(
#     aes(ymin=lower.ci.anom.rich,
#         ymax=upper.ci.anom.rich),
#     linetype=2, alpha=0.1, color="purple")+
#   theme(legend.position="none") +
#   facet_wrap(as.factor(RichAnom$system), scales = "free") +
#   geom_line()
# 
# # Density Deltas Distribution
# ggplot(DensDeltaClean, 
#        aes(x=deltaP)) + 
#   geom_histogram(binwidth=.5) +
#   facet_wrap(as.factor(DensDeltaClean$system), scales = "free")
# 
# ggplot(DensDeltaInf, 
#        aes(x=delta)) + 
#   geom_histogram(binwidth=.01) +
#   facet_wrap(as.factor(DensDeltaInf$system), scales = "free")

# Turnover
ggplot(data=TurnAnom,
       aes(x=year, y=total, color=system)) +
  geom_ribbon(
    aes(ymin=lower.ci.anom.total,
        ymax=upper.ci.anom.total),
    linetype=2, alpha=0.1, color="purple")+
  theme(legend.position="none") +
  theme_bw() + 
  facet_wrap(as.factor(TurnAnom$system), scales = "free") +
  geom_line()

ggplot(data=AppDis,
       aes(x=year, y=prop, color=system)) +
  #theme(legend.position="none") +
  facet_wrap(~direction) +
  geom_line(size = 2)

ggplot(data=total_turnover,
       aes(x=year, y=total, color=system)) +
  #theme(legend.position="none") +
  #facet_wrap(as.factor(AppDis$direction), scales = "free") +
  geom_line(size = 2)

ggplot(AllRate, aes(interval, distance, color = system)) + 
  facet_wrap(~system) + 
  geom_point() + 
  theme_bw() + 
  stat_smooth(method = "loess", se = F, size = 2)
