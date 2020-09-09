library(readxl)
library(tidyverse)
library(lubridate)

#Red Tide Switch
#Coded as months of interest: e.g., Jul-Oct
RTS <- c(7:10)

#NODC Codes to combine based on
#Difficult-to-ID taxa (DTI)
#AncSpp 8747020200
#EucSpp 8835390100

AncSpp <- c(8747020209, 8747020203, 8747020201, 8747020204, 8747020205, 8747020202, 8747020200)
EucSpp <- c(8835390101, 8835390102, 8835390111, 8835390108, 8835390109, 8835390104, 8835390105, 8835390100)

# import Excel files
NODCFull  <- read_excel("C:/Users/Gymnothorax/Google Drive/FIM_Data/gearCode20/TBM20AEZonesGridsRT.xlsx", sheet = "TBM20AEZones")
NODCCodes <- read_excel("C:/Users/Gymnothorax/Google Drive/FIM_Data/SpeciesCodes.xlsx", sheet = "SpeciesCodes")

# replace all DTI txa with string then replace that string with the appropriate NODC Code
# then group by reference and NODC Code to combine duplicates (other columns needed to show up after summarise)
# i.e., combines multiple reports of same NODC Code per reference (haul)
# covers modified codes above and multiple counts of taxa (as with LaRh)
NODC_NoDTI <- NODCFull %>%
  mutate(NODCCODE=replace(NODCCODE, NODCCODE %in% AncSpp, "AncSpp")) %>%
  mutate(NODCCODE=replace(NODCCODE, NODCCODE %in% EucSpp, "EucSpp")) %>%
  mutate(NODCCODE=replace(NODCCODE, NODCCODE == "AncSpp", 8747020200)) %>%
  mutate(NODCCODE=replace(NODCCODE, NODCCODE == "EucSpp", 8835390100)) %>%
  group_by(Reference, NODCCODE, Zone, Stratum, Grid, Sampling_Date) %>%
  summarise(Abundance = sum(Number))

# lookup NODC Codes against human readable names
NODC_to_HR <- merge(NODC_NoDTI, NODCCodes[,c(1:2)])

# select only RT months

# RT_Abund <- NODC_to_HR %>%
#   mutate(month = month(Sampling_Date), year = year(Sampling_Date)) %>%
#   mutate(RTLogic=(year == 2005)) %>%
#   filter(month %in% RTS)
RT_Abund <- NODC_to_HR %>%
  mutate(month = month(Sampling_Date), year = year(Sampling_Date)) %>%
  mutate(RTLogic = "Before") %>%
  mutate(RTLogic = replace(RTLogic, year == 2005, "During")) %>%
  mutate(RTLogic = replace(RTLogic, year > 2005, "After"))%>%
  filter(month %in% RTS)

# drop NODC Code here and then spread via scientific name
TrawlFull <- RT_Abund[,c(2:11)] %>%
  spread(Scientificname, Abundance) %>%
  replace(is.na(.),0)

# pull all environmental data out (basically not abundance data)
TrawlEnv <- TrawlFull %>%
  select(Reference:RTLogic) %>%
  as.data.frame()

# set zone as factor for rank abundance stuff
TrawlEnv$Zone <- as.factor(TrawlEnv$Zone)

# pull all abundance data out
TrawlAbun <- TrawlFull %>%
  subset(select = -c(Reference:RTLogic)) %>%
  as.data.frame()

#library(nlme)
#library(MASS)
library(vegan)
library(vegan3d)
#library(BiodiversityR)
library(goeveg)
library(scales)
library(ggrepel)

source('C:/Users/Gymnothorax/Box/Graduate/RGD/R Scripts/HaulWise.R')

##### SCHRAM - Rank Abundance Rebooted ####
Haul_PW <- HaulWise(TrawlFull, 1, "Zone", "RTLogic")

for (i in 1:length(Haul_PW)){
  temp_df <- as.data.frame(racurve(select(Haul_PW[[i]], -(Reference:RTLogic))))
  temp_df <- mutate(temp_df, Species = rownames(temp_df), .before = "abund")
  temp_df <- mutate(temp_df, Name = (rep(names(Haul_PW[i]), length(temp_df[[1]]))))
  temp_df <- separate(temp_df, Name, c("Zone","RTLogic"), sep="_")
  temp_df <- mutate(temp_df, FoO = (freq/nrow(Haul_PW[[i]])*100))
  temp_df <- mutate(temp_df, rank = seq.int(nrow(temp_df)))
  assign(paste(names(Haul_PW[i])), temp_df) 
}

Haul_RA <- rbind(A_After,A_Before,A_During, 
                 B_After,B_Before,B_During, 
                 C_After,C_Before,C_During, 
                 D_After,D_Before,D_During,
                 E_After,E_Before,E_During)

Haul_RA$RTLogic = factor(Haul_RA$RTLogic, levels = c("Before","During","After"))

rm(A_After,A_Before,A_During, 
   B_After,B_Before,B_During, 
   C_After,C_Before,C_During, 
   D_After,D_Before,D_During,
   E_After,E_Before,E_During)

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

##### END SCHRAM AMENDMENTS #####