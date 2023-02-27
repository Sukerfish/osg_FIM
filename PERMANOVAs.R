# PERMANOVAs for each system/season combination


#### data input ######
library(tidyverse)
library(vegan)

load('TidyGearCode20.Rdata')

#get the biological data and associated site chars
CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, BottomVegCover, systemZone, Scientificname, N2))
#filter(season == "summer" | season == "winter")
#filter(season == "winter")
#filter(season == "summer")

#get water temp from hydrology dataset
WaterTemp <- HydroList %>%
  subset(select = c(Reference, Temperature))

#collect everything by Reference and spread to long form
SXS_full <- CleanHauls %>%
  mutate(N2 = N2^.25) %>% #fourth-root transform
  group_by(Reference) %>%
  filter(sum(N2)>0) %>% #remove all References with 0 taxa found
  spread(Scientificname,N2) %>%
  ungroup() %>%
  #subset(select = -c(Reference, season, systemZone)) %>%
  replace(is.na(.), 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  group_by(season, system, seasonYear) %>%
  mutate(seasonYear = as.factor(seasonYear)) %>%
  left_join(WaterTemp) %>%
  filter(!is.na(Temperature)) %>% #only excludes 7 sampling events from above pool
  unite("systemSeason",
        system:season,
        sep = "_") %>%
  ungroup()
#summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))


##### setup for PERMANOVA loop ####
systemSeason_list <- SXS_full %>%
  select(systemSeason) %>%
  distinct()

PERMSforAll <- list()
for(i in systemSeason_list$systemSeason){
  print(i) #watch progress through list
  
  df <- SXS_full %>% #filter out systemSeason of interest
    filter(systemSeason %in% i)

  df_spe <- df %>% #pull out taxa only
    subset(select = -c(systemSeason, seasonYear, Reference, systemZone, BottomVegCover, Temperature)) %>%
    select(which(!colSums(., na.rm=TRUE) %in% 0))

  df_env <- df %>% #pull out environmental variables
    subset(select = c(systemSeason, seasonYear, BottomVegCover, Temperature)) %>%
    mutate(contYear = as.numeric(as.character(seasonYear)))

  bf <- (vegdist(df_spe))^0.5 #Bray-Curtis w/ sqrt to prevent negative eigenvalues
  
  PERM = adonis2(bf ~ contYear * Temperature + BottomVegCover,
                              data = df_env,
                              #add = "lingoes",
                              parallel = 6,
                              #method="bray", 
                              permutations = 999)
  
  PERMSforAll[[i]] <- PERM
}

PERMSforAllDF <- bind_rows(PERMSforAll, .id = "systemSeason") %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_") %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))

write.csv(PERMSforAllDF, './Outputs/PERMSforALLDF.csv')