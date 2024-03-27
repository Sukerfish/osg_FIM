# PERMANOVAs for each system/season combination


#### data input ######
library(tidyverse)
library(vegan)

load('TidyGearCode20.Rdata')
#load('SXS_filtered.Rdata')
load("SXS_filtered_fars.Rdata")

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


#build month dataframe
monthly <- HydroList %>%
  select(Reference, month)

#rejoin with environmental data
SXS_run_env <- SXS_filtered_env %>%
  #left_join(SXS_filtered_env) %>%
  left_join(monthly) %>%
  group_by(systemSeason) %>%
  add_count(name = "n_hauls") %>% #count number of hauls per systemSeason
  ungroup() %>%
  group_by(systemSeason, seasonYear) %>%
  mutate(avg_temp = mean(Temperature, na.rm = TRUE),
         n_temp = n(),
         upper = quantile(Temperature, 0.9),
         lower = quantile(Temperature, 0.1),
         sd_ann = sd(Temperature)) %>%
  ungroup() %>%
  group_by(systemSeason) %>%
  mutate(avg_ltm = mean(avg_temp),
         sd_ltm = sd(avg_temp, na.rm = TRUE),
         upper_sd = sd(upper),
         lower_sd = sd(lower)) %>%
  mutate(anom_temp = avg_temp - avg_ltm,
         se_temp = sd_ltm/sqrt(n_temp),
         lower.ci.anom.temp = 0 - (1.96 * se_temp),
         upper.ci.anom.temp = 0 + (1.96 * se_temp)) %>%
  #zscoreing
  mutate(bvc_ltm  = mean(BottomVegCover, na.rm = TRUE),
         bvc_sd   = sd(BottomVegCover, na.rm = TRUE),
         temp_ltm = mean(Temperature, na.rm = TRUE),
         temp_sd  = sd(Temperature, na.rm = TRUE),
         year_ltm = mean(as.numeric(as.character(seasonYear)), na.rm = TRUE),
         year_sd  = sd(as.numeric(as.character(seasonYear)), na.rm = TRUE),
         bvc_Z    = ((BottomVegCover - bvc_ltm)/bvc_sd),
         temp_Z   = ((Temperature - temp_ltm)/temp_sd),
         year_Z    = ((as.numeric(as.character(seasonYear)) - year_ltm)/year_sd),
         sd_t_Z   = sd(temp_Z)) %>%
  #monthly
  ungroup() %>%
  group_by(systemSeason, seasonYear, month) %>%
  mutate(avg_temp_mon = mean(Temperature, na.rm = TRUE),
         n_temp_mon = n(),
         upper_mon = quantile(Temperature, 0.9),
         lower_mon = quantile(Temperature, 0.1),
         upper_Z  = quantile(temp_Z, 0.9),
         lower_Z  = quantile(temp_Z, 0.1),
         sd_mon = sd(Temperature),
         se_mon = sd_mon/sqrt(n_temp_mon))

##### setup for PERMANOVA loop ####
systemSeason_list <- SXS_full %>%
  select(systemSeason) %>%
  distinct()

PERMSforAll <- list()
for(i in systemSeason_list$systemSeason){
  print(i) #watch progress through list
  
  df <- SXS_filtered %>% #filter out systemSeason of interest
    filter(systemSeason %in% i)
  
  # df_spe <- df %>% #pull out taxa only
  #   subset(select = -c(systemSeason, seasonYear, Reference, systemZone, BottomVegCover, Temperature)) %>%
  #   select(which(!colSums(., na.rm=TRUE) %in% 0)) #select only taxa present in this systemSeason
  # 
  # df_pa <- df_spe
  # df_pa[df_pa > 0] <- 1 #convert to pa
  # 
  # spp <- length(df_spe)
  # spx <- nrow(df_spe)
  # 
  # df_pa_filtered <- df_pa %>%
  #   select_if(colSums(.)>(0.05*spx))
  # 
  # df_filtered <- df %>%
  #   select(c(Reference, seasonYear, all_of(colnames(df_pa_filtered)))) %>%
  #   rowwise() %>%
  #   mutate(N = sum(across(!c(Reference:seasonYear)))) %>%
  #   ungroup() %>%
  #   filter(N > 0) %>%
  #   select(!c(N))
  
  df_spe_filtered <- SXS_filtered %>%
    filter(Reference %in% df$Reference) %>%
    subset(select = -c(systemSeason, seasonYear, Reference, systemZone, BottomVegCover, Temperature)) %>%
    select(which(!colSums(., na.rm=TRUE) %in% 0)) #select only taxa present in this systemSeason
  
  df_env <- data.frame(SXS_run_env %>% #pull out environmental variables
                         filter(Reference %in% df$Reference) %>%
                         subset(select = c(systemSeason, seasonYear, bvc_Z, temp_Z)) %>%
                         mutate(contYear = as.numeric(as.character(seasonYear))))
  
  bf <- (vegdist(df_spe_filtered))^0.5 #Bray-Curtis w/ sqrt to reduce negative eigenvalues
  
  PERM = adonis2(bf ~ seasonYear,
                              data = df_env,
                #strata = df_env$seasonYear,
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

#save(PERMSforAllDF, file = './Outputs/PERMSforALLDF_factor_fars.RData')
