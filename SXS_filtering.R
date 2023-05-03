# #### Import and Filter Data #####

library(tidyverse)

load('TidyGearCode20.Rdata')

#get the biological data and associated site chars
CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, BottomVegCover, systemZone, Scientificname, N2))

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

##### set up loop ####
systemSeason_list <- SXS_full %>%
  select(systemSeason) %>%
  distinct()

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

  spp <- length(df_spe) #how many taxa total
  spx <- nrow(df_spe) #how many samples present

  df_pa_filtered <- df_pa %>%
    select_if(colSums(.)>(0.05*spx)) #keep only columns (taxa) that have a total abundance of >5% (i.e., not rare)

  df_filtered <- df %>%
    select(c(Reference, seasonYear, all_of(colnames(df_pa_filtered)))) %>% #start again with only common taxa
    rowwise() %>%
    mutate(n = sum(across(!c(Reference:seasonYear)))) %>% #determine if any references had ONLY rare taxa
    ungroup() %>%
    filter(n > 0) %>% #filter out references that were ONLY rare taxa
    select(Reference)

SXS_filteredList[[i]] <- df_filtered #export per systemSeason

}

SXS_filtered <- bind_rows(SXS_filteredList) %>%
  left_join(SXS_full) #join back into main SXS data but with only filtered refs

#pull out relevant subsections
SXS_filtered_env <- SXS_filtered %>%
  select(c(Reference, systemSeason, seasonYear, systemZone, BottomVegCover, Temperature))

SXS_filtered_spp <- SXS_filtered %>%
  select(!c(Reference, systemSeason, seasonYear, systemZone, BottomVegCover, Temperature))

#save(SXS_filtered, SXS_filtered_env, SXS_filtered_spp, file = "SXS_filtered.Rdata")