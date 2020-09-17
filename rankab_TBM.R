library(tidyverse)
library(lubridate)
library(RODBC)

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

# clean up biology and merge with selected references --------

CleanBio20 <- BioNumFull %>%
  select(Reference, Species_record_id, Splittype, Splitlevel, Cells, 
         NODCCODE, Number, FHC) %>%
  filter(FHC != "D" &
           #remove unidentified species
           !NODCCODE %in% c("9999000000") &
           #remove freshwater prawns
           !NODCCODE %in% c("6179110200", "6179110201")) %>%
  collect() %>%
  filter(#remove all turtles, by NODCCODEs that start with 9002
    !str_detect(NODCCODE, "^9002"),
    #remove all grass shrimp, by NODCCODEs that start with 617911
    !str_detect(NODCCODE, "^617911")) %>%
  #subset to references of interest
  inner_join(RefsList, by = "Reference") %>%
  #account for splitter
  mutate(Count = case_when(!is.na(as.numeric(Splittype)) ~ Number*(as.numeric(Splittype)^as.numeric(Splitlevel)),
                           TRUE ~ as.numeric(Number))) %>%
  select(-Splittype, -Splitlevel, -Species_record_id, -Cells, -FHC) %>%
  group_by(Reference, NODCCODE) %>%
  mutate(N = sum(Count)) %>%
  select(-Count, -Number) %>%
  distinct() %>%
  ungroup() %>%
  #select only a few variables and expand sampling date
  select(Reference, NODCCODE, Sampling_Date, Stratum, Zone,
         Grid, BottomVegCover,BycatchQuantity, Bank, ShoreDistance, N) %>%
  mutate(month = month(Sampling_Date), year = year(Sampling_Date)) %>%
  select(-Sampling_Date)

# -----------

#Red Tide Switch
#Coded as months of interest: e.g., Jul-Oct
RTS <- c(7:10)
StartYear <- 1998
EndYear <- 2017

#NODC Codes to combine based on
#Difficult-to-ID taxa (DTI)
#AncSpp 8747020200
#EucSpp 8835390100

AncSpp <- c(8747020209, 8747020203, 8747020201, 8747020204, 8747020205, 8747020202, 8747020200)
EucSpp <- c(8835390101, 8835390102, 8835390111, 8835390108, 8835390109, 8835390104, 8835390105, 8835390100)

# replace all DTI txa with string then replace that string with the appropriate NODC Code
# then group by reference and NODC Code to combine duplicates
# i.e., combines multiple reports of same NODC Code per reference (haul)
# covers modified codes above and multiple counts of taxa (as with LaRh)
NODC_NoDTI <- CleanBio20 %>%
  mutate(NODCCODE=replace(NODCCODE, NODCCODE %in% AncSpp, "AncSpp")) %>%
  mutate(NODCCODE=replace(NODCCODE, NODCCODE %in% EucSpp, "EucSpp")) %>%
  mutate(NODCCODE=replace(NODCCODE, NODCCODE == "AncSpp", 8747020200)) %>%
  mutate(NODCCODE=replace(NODCCODE, NODCCODE == "EucSpp", 8835390100)) %>%
  group_by(Reference, NODCCODE) %>%
  mutate(N2 = sum(N)) %>%
  select(-N) %>%
  distinct() %>%
  ungroup()

# lookup NODC Codes against human readable names
CleanHRBio <- inner_join(NODC_NoDTI, SpeciesList, by = "NODCCODE") %>%
  select(Reference, Scientificname, month, year, Stratum, Zone,
         Grid, BottomVegCover,BycatchQuantity, #Bank, 
         ShoreDistance, N2) %>%
  mutate(#normalize Zone strings
    Zone = str_trim(str_to_upper(Zone))) %>%
  mutate(#normalize Stratum strings
    Stratum = str_trim(str_to_upper(Stratum)))

# select only RT months
RT_Abund <- CleanHRBio %>%
  mutate(RTLogic = "Before") %>%
  mutate(RTLogic = replace(RTLogic, year == 2005, "During")) %>%
  mutate(RTLogic = replace(RTLogic, year > 2005, "After")) %>%
  filter(month %in% RTS) %>%
  filter(#remove early years
    year >= StartYear)

# spread via scientific name
HaulFull <- RT_Abund %>%
  spread(Scientificname, N2) %>%
  replace(is.na(.),0)

# pull all environmental data out
HaulEnv <- HaulFull %>%
  select(Reference:RTLogic) %>%
  as.data.frame()

# set zone as factor for rank abundance stuff
HaulEnv$Zone <- as.factor(HaulEnv$Zone)

# pull all abundance data out
HaulAbun <- HaulFull %>%
  subset(select = -c(Reference:RTLogic)) %>%
  as.data.frame()

# summary information about the samples
HaulCount <- HaulFull %>%
  count(year, Zone)
HaulMin <- min(HaulCount$n)

library(vegan)
library(vegan3d)
library(goeveg)
library(scales)
library(ggrepel)
#library(devtools)

#source_url("https://github.com/Sukerfish/osg_FIM/blob/master/HaulWise.R?raw=TRUE")

##### Yearly Resampled Species Richness ####

# Resample using the minimum samples in the paired values
HaulSub <- HaulFull %>%
  group_by(year, Zone) %>%
  sample_n(HaulMin, replace = FALSE) %>%
  ungroup()

Haul_Rich <- data.frame("year" = double(),
                        "A" = double(),
                        "B" = double(),
                        "C" = double(),
                        "D" = double(),
                        "E" = double())

for (i in 1:20){
  temp_df <- HaulSub %>%
    filter(year == i+(StartYear-1))
  Haul_Rich[i,1]  <- i+(StartYear-1)
  Haul_Rich[i,-1] <- specnumber(temp_df, temp_df$Zone)
}

Haul_Rich <- Haul_Rich %>%
  gather('A','B','C','D','E', key = "Zone", value = "Richness")

####### OLD METHODS #######
# Haul_PW <- HaulWise(HaulFull, 1, "year", "Zone")
# Haul_PW_list <- names(Haul_PW)
# 
# for (i in 1:length(Haul_PW)){
#   temp_rs <- select(# grab only the abundances for each
#     Haul_PW[[i]], -(Reference:RTLogic))
#   # resample using the minimum sample size from the above pairs
#   temp_rs <- temp_rs[sample(nrow(temp_rs), HaulMin), ]
#   # calculate rankabundance data
#   temp_df <- as.data.frame(racurve(temp_rs))
#   temp_df <- mutate(# get species names from rows
#     temp_df, Species = rownames(temp_df), .before = "abund")
#   
#   temp_df <- mutate(temp_df, Name = (rep(names(Haul_PW[i]), length(temp_df[[1]]))))
#   temp_df <- separate(temp_df, Name, c("year","Zone"), sep="_")
#   temp_df <- mutate(temp_df, FoO = (freq/nrow(temp_rs)*100))
#   temp_df <- mutate(temp_df, rank = seq.int(nrow(temp_df)))
#   assign(paste(names(Haul_PW[i])), temp_df) 
# }
# 
# Haul_RA <- rbind(A_After,A_Before,A_During, 
#                  B_After,B_Before,B_During, 
#                  C_After,C_Before,C_During, 
#                  D_After,D_Before,D_During,
#                  E_After,E_Before,E_During)
# 
# Haul_RA$RTLogic = factor(Haul_RA$RTLogic, levels = c("Before","During","After"))
# 
# rm(A_After,A_Before,A_During, 
#    B_After,B_Before,B_During, 
#    C_After,C_Before,C_During, 
#    D_After,D_Before,D_During,
#    E_After,E_Before,E_During)
# 
# rm(Haul_PW)
# 
# 
# Haul_PW <- HaulWise(HaulFull, 1, "Zone", "RTLogic")
# 
# for (i in 1:length(Haul_PW)){
#   temp_df <- as.data.frame(racurve(select(Haul_PW[[i]], -(Reference:RTLogic))))
#   temp_df <- mutate(temp_df, Species = rownames(temp_df), .before = "abund")
#   temp_df <- mutate(temp_df, Name = (rep(names(Haul_PW[i]), length(temp_df[[1]]))))
#   temp_df <- separate(temp_df, Name, c("Zone","RTLogic"), sep="_")
#   temp_df <- mutate(temp_df, FoO = (freq/nrow(Haul_PW[[i]])*100))
#   temp_df <- mutate(temp_df, rank = seq.int(nrow(temp_df)))
#   assign(paste(names(Haul_PW[i])), temp_df) 
# }
# 
# Haul_RA <- rbind(A_After,A_Before,A_During, 
#                  B_After,B_Before,B_During, 
#                  C_After,C_Before,C_During, 
#                  D_After,D_Before,D_During,
#                  E_After,E_Before,E_During)
# 
# Haul_RA$RTLogic = factor(Haul_RA$RTLogic, levels = c("Before","During","After"))
# 
# rm(A_After,A_Before,A_During, 
#    B_After,B_Before,B_During, 
#    C_After,C_Before,C_During, 
#    D_After,D_Before,D_During,
#    E_After,E_Before,E_During)
# 
# rm(Haul_PW)
####### Plots ######

#Beta diversity over time plot
ggplot(Haul_Rich, aes(x=year))+
  geom_line(aes(y=Richness, group = 1))+
  theme_bw()+
  labs(title = "Richness Over Time")+
  theme(axis.title = element_text(face = "bold", size = "12"), 
        axis.text = element_text(size = "12"),
        strip.text = element_text(face = "bold", size = "14"))+
  scale_x_continuous(name = "Year", 
                     limits = c(StartYear,EndYear),
                     breaks = seq(StartYear,EndYear, 2))+
  scale_y_continuous(name = "Number of Taxa",
                     limits = c(0,70))+
  theme(legend.position="bottom")+
  facet_grid(Zone ~ .)

#Ranked Abundance Plot
ggplot(Haul_RA, aes(rank, abund))+
  geom_point()+
  geom_line()+
  theme_bw()+
  labs(title = "Ranked abundance curves")+
  theme(axis.title = element_text(face = "bold", size = "12"), 
        axis.text = element_text(size = "12"),
        strip.text = element_text(face = "bold", size = "14"))+
  scale_x_continuous(name = "", limits = c(0,110), breaks = seq(0,108,12))+
  scale_y_continuous(name = "Log10(Abundance)",
                     trans = log10_trans(),
                     #limits = c(10^-3, 10^5),
                     #breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  geom_text_repel(aes(label=ifelse(rank <= 5, as.character(Species), "")),
                  force = 1,
                  size = 3.5,
                  nudge_x = 95,
                  direction = "y",
                  segment.color = "grey"
  )+
  facet_grid(Zone~RTLogic)

#Ranked Frequency of Occurence Plot
ggplot(Haul_RA, aes(rank, FoO))+
  geom_point()+
  geom_line()+
  theme_bw()+
  labs(title = "Ranked frequency of occurence curves")+
  theme(axis.title = element_text(face = "bold", size = "12"), 
        axis.text = element_text(size = "12"),
        strip.text = element_text(face = "bold", size = "14"))+
  scale_x_continuous(name = "", limits = c(0,110), breaks = seq(0,108,12))+
  scale_y_continuous(name = "Frequency of occurrence")+
  geom_text_repel(aes(label=ifelse(rank <= 5, as.character(Species), "")),
                  force = 1,
                  size = 3.5,
                  nudge_x = 30,
                  direction = "y",
                  segment.color = "grey")+
  facet_grid(Zone~RTLogic)