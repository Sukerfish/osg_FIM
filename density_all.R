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
  filter(str_detect(Reference,"^APM"))

# Establish Zone filter
ZoneFilter = c("A"
               ,"B"
               #,"C"
               #,"D"
               #,"E"
)

# clean up biology and merge with selected references
CleanBio   <- osg_CleanBio(BioNumFull, RefsList)

# merge difficult to ID taxa and convert NODC code to human readable
CleanHRBio <- osg_ComBio(CleanBio, SpeciesList)

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
  mutate(season = "Winter") %>%
  mutate(season = replace(season, month %in% c(4:7), "Spring")) %>%
  mutate(RTLogic = "Before") %>%
  mutate(RTLogic = replace(RTLogic, year == 2005, "During")) %>%
  mutate(RTLogic = replace(RTLogic, year > 2005, "After")) %>%
  filter(month %in% RTS) %>%
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

##### Yearly Resampled Species density ####

# Resample using the minimum samples in the paired values
HaulSub <- HaulFullCleanTBM %>%
  group_by(year) %>%
  # resample point
  sample_n(HaulMin, replace = FALSE) %>%
  ungroup()

HaulSub_Abun <- HaulSub %>%
  subset(select = -c(Stratum:RTLogic)) %>%
  subset(select = -c(total:CPUE)) %>%
  subset(select = -c(Reference:month))

HaulSub_Dens <- HaulSub_Abun %>%
  mutate(across(!c(year), ~{.x/effort}))

# # initialize
# testing <- setNames(data.frame(matrix(ncol = ncol(HaulSub_Abun),
#                                       nrow = 0)),
#                     c(colnames(HaulSub_Abun)))
# for (i in 1:((EndYear-StartYear)+1)){
#   temp_df <- HaulSub_Abun %>%
#     # filter by year and remove extra columns
#     filter(year == i+(StartYear-1))
#     #subset(select = CPUE)
#   # assign year to row
#   testing[i,1]  <- i+(StartYear-1)
#   # calculate richness per haul and average it before placement in matrix
#   for (j in 2:ncol(HaulSub_Abun)){
#     testing[i,j] <- mean(temp_df[,j], na.rm = TRUE)
#   }
# }

testtt <- HaulSub_Dens %>%
  group_by(year) %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(!year, names_to = "species", values_to = "CPUE")



# calculate average density per haul using the subsampled data
# initialize
Haul_AvgDensity <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("year", "avg_density"))
for (i in 1:((EndYear-StartYear)+1)){
  temp_df <- HaulSub %>%
    # filter by year and remove extra columns
    filter(year == i+(StartYear-1)) %>%
    subset(select = CPUE)
    # assign year to row
    Haul_AvgDensity[i,1]  <- i+(StartYear-1)
    # calculate richness per haul and average it before placement in matrix
    Haul_AvgDensity[i,2] <- mean(temp_df$CPUE)
  }

# summary stats and plotting processing
Haul_tsDensity <- Haul_AvgDensity %>%
  # calculate as anomalies
  mutate(anom.Dens = avg_density - mean(avg_density, na.rm = TRUE))

DensityAnom_Metrics <- Haul_tsDensity %>%
  # calculate mean/CIs
  summarise(mean.Dens = mean(avg_density, na.rm = TRUE),
            sd.anom.Dens = sd(anom.Dens, na.rm = TRUE),
            n.Dens = n()) %>%
  mutate(se.anom.Dens = sd.anom.Dens / sqrt(n.Dens),
         lower.ci.anom.Dens = 0 - (1.96 * se.anom.Dens),
         upper.ci.anom.Dens = 0 + (1.96 * se.anom.Dens))

####### Plots ######
# Richness over time plot
ggplot(Haul_tsDensity, aes(x=year))+
  geom_line(aes(y=anom.Dens))+
  geom_hline(aes(yintercept = 0, color = "darkred"))+
  geom_ribbon(data=merge(DensityAnom_Metrics, Haul_tsDensity),
              aes(ymin=lower.ci.anom.Dens, 
                  ymax=upper.ci.anom.Dens), 
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


ggplot(data=testtt,
       aes(x=year, y=CPUE, color=species)) +
  theme(legend.position="none") +
  facet_wrap(as.factor(testtt$species), scales = "free") +
  geom_line()


ggplot(testtt, aes(x=year))+
  geom_line(aes(y=CPUE))+
  #geom_hline(aes(yintercept = 0, color = "darkred"))+
  #geom_ribbon(data=merge(DensityAnom_Metrics, Haul_tsDensity),
              #aes(ymin=lower.ci.anom.Dens, 
              #    ymax=upper.ci.anom.Dens), 
             # linetype=2, alpha=0.1, color="purple")+
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