# Bray-Curtis similarity pairwise among all years with respect to first year of data
# Using average annual density of each taxa within each estuary and season
# 
# Model change in abundance through time for each taxa
# 
# Calculate coefficient of the slope of each
# 
# Potential coefficient distributions
#### data input ######

library(tidyverse)
library(vegan)

load('TidyGearCode20.Rdata')

CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, BottomVegCover, systemZone, Scientificname, N2))
#filter(season == "summer" | season == "winter")
#filter(season == "winter")
#filter(season == "summer")

WaterTemp <- HydroList %>%
  subset(select = c(Reference, Temperature))

###### PERMANOVAs pairwise ######
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

systemSeason_list <- SXS_full %>%
  select(systemSeason) %>%
  distinct()

## filter out hauls that were solely rare (<=5% presence in total samples)
SXS_filteredList <- list()
for(i in systemSeason_list$systemSeason){
  print(i) #watch progress through list
  
  df <- SXS_full %>% #filter out systemSeason of interest
    filter(systemSeason %in% i)
  
  df_spe <- df %>% #pull out taxa only
    subset(select = -c(systemSeason, seasonYear, Reference, systemZone, BottomVegCover, Temperature)) %>%
    select(which(!colSums(., na.rm=TRUE) %in% 0)) #select only taxa present in this systemSeason
  
  df_pa <- df_spe
  df_pa[df_pa > 0] <- 1 #convert to pa
  
  spp <- length(df_spe)
  spx <- nrow(df_spe)
  
  df_pa_filtered <- df_pa %>%
    select_if(colSums(.)>(0.05*spx))
  
  df_filtered <- df %>%
    select(c(Reference, seasonYear, all_of(colnames(df_pa_filtered)))) %>%
    rowwise() %>%
    mutate(n = sum(across(!c(Reference:seasonYear)))) %>%
    ungroup() %>%
    filter(n > 0) %>%
    select(Reference)
  
  SXS_filteredList[[i]] <- df_filtered
  
}

SXS_filtered2 <- bind_rows(SXS_filteredList) %>%
  left_join(SXS_full) %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_") %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH"))) %>%
  mutate(seasonYear = as.factor(seasonYear)) %>%
  filter(system == "CK") %>%
  mutate(event = ifelse(as.numeric(as.character(seasonYear)) < 2011, "pre", "post"))

CK_summer_spe <- SXS_filtered2 %>%
  filter(season == "winter") %>%
  subset(select = -c(season, system, seasonYear, Reference,systemZone, BottomVegCover,season, Temperature, event)) %>%
  select(which(!colSums(., na.rm=TRUE) %in% 0))

CK_summer_env <- SXS_filtered2 %>%
  filter(season == "winter") %>%
  subset(select = c(system,seasonYear,BottomVegCover, Temperature, event)) %>%
  mutate(contYear = as.numeric(as.character(seasonYear))) %>%
  mutate(smallYear = contYear - min(contYear))

# bray curtis distance... square rooted - to make all eigenvalues positive per 
CK_summer_bray = vegdist(CK_summer_spe)
CK_summer_bray2 = CK_summer_bray^.5

# homdispSyWinter = permutest(betadisper(SXSf_winter_bray2, SXSf_winter_env$system, type = "centroid"), 
#                             pairwise = TRUE, permutations = 999, parallel = 6)
# 
# homdispYrWinter = permutest(betadisper(SXSf_winter_bray2, SXSf_winter_env$seasonYear, type = "centroid"),
#                             pairwise = TRUE, permutations = 999, parallel = 6)


postTest = adonis2(CK_summer_bray2 ~ event,
                           data = CK_summer_env,
                           #add = "lingoes",
                           parallel = 6,
                           #method="bray",
                           permutations=999)
