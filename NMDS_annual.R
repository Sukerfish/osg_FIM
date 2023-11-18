##### NMDS annual abundance data #####

library(tidyverse)
library(vegan)
library(pixiedust)

load('SXS_filtered.RData')

monthly <- HydroList %>%
  select(Reference, month)

#### winter first ####
test_winter <- SXS_filtered %>%
  left_join(monthly) %>%
  filter(systemSeason %in% c("AP_winter", "CK_winter", "CH_winter", "TB_winter")) %>%
  group_by(systemSeason, seasonYear) %>%
  select(!c(Reference, systemZone)) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

test_env_winter <- test_winter %>%
  select(systemSeason, seasonYear, month, BottomVegCover, Temperature)

test_spp_winter <- test_winter %>%
  ungroup() %>%
  select(!c(systemSeason, seasonYear, month, BottomVegCover, Temperature)) %>%
  select(which(!colSums(., na.rm=TRUE) %in% 0))

nmds_ann_winter = metaMDS(test_spp_winter,
                      distance = "bray",
                      k=2,
                      maxit = 50,
                      trymax = 50,
                      wascores = TRUE)

env_winter = envfit(nmds_ann_winter, test_env_winter, permutations = 999, na.rm = TRUE)

plot(nmds_ann_winter, type = "n", las = 1, main = "NMDS winter")
points(nmds_ann_winter, display = "sites")
points(nmds_ann_winter, display = "species", col = "red", pch = 3)
ordiellipse(nmds_ann_winter,
            groups = test_env_winter$systemSeason,
            kind = "se",
            conf = 0.95,
            display = "sites",
            label=T)

plot(nmds_ann_winter, main = "NMDS winter")
ordihull(nmds_ann_winter,groups=test_env_winter$systemSeason,draw="polygon",col="grey90",label=F)
#plot(SXSRaw_sppfit, p.max = 0.001, col = "black", cex = 0.7)

#### summer next ####
test_summer <- SXS_filtered %>%
  left_join(monthly) %>%
  filter(systemSeason %in% c("AP_summer", "CK_summer", "CH_summer", "TB_summer")) %>%
  group_by(systemSeason, seasonYear) %>%
  select(!c(Reference, systemZone)) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

test_env_summer <- test_summer %>%
  select(systemSeason, seasonYear, month, BottomVegCover, Temperature)

test_spp_summer <- test_summer %>%
  ungroup() %>%
  select(!c(systemSeason, seasonYear, month, BottomVegCover, Temperature)) %>%
  select(which(!colSums(., na.rm=TRUE) %in% 0))

nmds_ann_summer = metaMDS(test_spp_summer,
                          distance = "bray",
                          k=2,
                          maxit = 50,
                          trymax = 50,
                          wascores = TRUE)

env_summer = envfit(nmds_ann_summer, test_env_summer, permutations = 999, na.rm = TRUE)

plot(nmds_ann_summer, type = "n", las = 1, main = "NMDS summer")
points(nmds_ann_summer, display = "sites")
points(nmds_ann_summer, display = "species", col = "red", pch = 3)
ordiellipse(nmds_ann_summer,
            groups = test_env_summer$systemSeason,
            kind = "se",
            conf = 0.95,
            display = "sites",
            label=T)

plot(nmds_ann_summer, main = "NMDS summer")
ordihull(nmds_ann_summer,groups=test_env_summer$systemSeason,draw="polygon",col="grey90",label=F)
#plot(SXSRaw_sppfit, p.max = 0.001, col = "black", cex = 0.7)