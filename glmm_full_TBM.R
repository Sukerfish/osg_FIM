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
HaulFull$Zone    <- as.factor(HaulFull$Zone)
HaulFull$RTLogic <- as.factor(HaulFull$RTLogic)
#HaulFull$RTLogic <- ordered(HaulFull$RTLogic, levels = c("Before", "During", "After"))
HaulFull$Stratum <- as.factor(HaulFull$Stratum)

# pull all environmental data out
HaulEnv <- HaulFull %>%
  select(Reference:RTLogic) %>%
  as.data.frame()
# 
# # pull all abundance data out
# HaulAbun <- HaulFull %>%
#   subset(select = -c(Reference:RTLogic)) %>%
#   as.data.frame()
# 
# # summary information about the samples
# HaulCount <- HaulFull %>%
#   count(year, Zone)
# HaulMin <- min(HaulCount$n)
# ZoneCount <- length(ZoneFilter)

BottomVegMetrics <- HaulFull %>%
  #select(-year) %>%
  group_by(Zone) %>%
  #drop_na() %>%
  # calculate mean/CIs
  summarise(mean.Cov = mean(BottomVegCover),
            sd.Cov = sd(BottomVegCover),
            n.Cov = n(),
            #max.Cov = max(BottomVegCover),
            #min.Cov = min(BottomVegCover),
            ) %>%
  mutate(se.Cov = sd.Cov / sqrt(n.Cov),
         lower.ci.Cov = mean.Cov - (1.96 * se.Cov),
         upper.ci.Cov = mean.Cov + (1.96 * se.Cov))

library(vegan)
#library(vegan3d)
#library(goeveg)
#library(scales)
#library(ggrepel)
library(MASS)
library(nlme)

##### Mixed Effects model ####

# # Resample using the minimum samples in the paired values
# HaulSub <- HaulFull %>%
#   group_by(year, Zone) %>%
#   # resample point
#   sample_n(HaulMin, replace = FALSE) %>%
#   ungroup()

HaulRich <- HaulFull %>%
  subset(select = -c(month:RTLogic))

HaulRich$Richness <- specnumber(HaulRich, HaulRich$Reference)

HaulRich <- HaulRich %>%
  subset(select = c(Reference, Richness)) %>%
  merge(HaulEnv)

HaulAbun <- RT_Abund %>%
  group_by(Reference) %>%
  summarise(totalAb = sum(N2)) %>%
  merge(HaulEnv)
  
modelAbun <- glmmPQL(totalAb~RTLogic * Zone,
                 random=~1|Grid,
                 family = "poisson",
                 data = HaulAbun,
                 corr=corARMA(form=~1|Grid/year,p=1),
                 )
summary(modelAbun)

modelRich <- glmmPQL(Richness~RTLogic * Zone,
                 random=~1|Grid,
                 family = "poisson",
                 data = HaulRich,
                 corr=corARMA(form=~1|Grid/year,p=1),
                 )
summary(modelRich)

modelAbunSAV <- glmmPQL(totalAb~RTLogic * Zone + BottomVegCover + SAVsq,
                     random=~1|Grid,
                     family = "poisson",
                     data = HaulAbun,
                     corr=corARMA(form=~1|Grid/year,p=1),
                     )
summary(modelAbunSAV)

modelRichSAV <- glmmPQL(Richness~RTLogic * Zone + BottomVegCover + SAVsq,
                     random=~1|Grid,
                     family = "poisson",
                     data = HaulRich,
                     corr=corARMA(form=~1|Grid/year,p=1),
)
summary(modelRichSAV)

# tickdata = read.table("~/Elston2001_tickdata.txt",header=TRUE,
#                       colClasses=c("factor","numeric","factor","numeric","factor","factor"))
# tickdata <- transform(tickdata,cHEIGHT=HEIGHT-mean(HEIGHT))
# 
# tmod_PQL <- glmmPQL(TICKS~cHEIGHT+YEAR,
#                     random=~1|LOCATION/BROOD,
#                     family="poisson",data=tickdata,
#                     verbose=FALSE)
# 
# ggplot(tickdata,aes(x=HEIGHT,y=1+TICKS,colour=YEAR))+
#   stat_sum(aes(size=..n..),alpha=0.7)+
#   scale_y_log10()+
#   scale_size_continuous(breaks=c(2,6,10),range=c(2,7))+
#   geom_smooth(method="glm",method.args=list(family=quasipoisson))

ggplot(test,aes(x=RTLogic,y=1+totalAb,colour=Zone))+
  stat_sum(aes(size=..n..),alpha=0.7)+
  scale_y_log10()+
  scale_size_continuous(breaks=c(5,10,15),range=c(2,7))+
  geom_smooth(method="glm",method.args=list(family=quasipoisson))