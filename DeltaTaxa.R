# Heatmap of zscored abundance data
#
#
#### data input ######

library(viridis)
library(tidyverse)
library(ggplot2)

load('TidyGearCode20.Rdata')

#### richness ####
CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, systemZone, Scientificname, N2))

# summerHauls <- CleanHauls %>%
#   filter(season == "summer")
# 
# winterHauls <- CleanHauls %>%
#   filter(season == "winter")
# 
# rawAbs <- CleanHauls %>%
#   subset(select = c(Reference, system, season, seasonYear, Scientificname, N2))
# 
# means <- CleanHauls %>%
#   #mutate(N2 = N2^.25) %>% #fourth-root transform
#   group_by(system, season, Scientificname) %>%
#   summarise(mean = mean(N2))
# 
# stdev <- CleanHauls %>%
#   group_by(system, season, Scientificname) %>%
#   summarise(stdev = sd(N2))
# 
# centered <- rawAbs %>%
#   group_by(system, season, seasonYear, Scientificname) %>%
#   summarise(avg = mean(N2)) %>%
#   left_join(means) %>%
#   left_join(stdev) %>%
#   mutate(zscore = ((avg - mean)/stdev)) %>%
#   filter(Scientificname != "No fish") %>%
#   filter(Scientificname == "Opsanus beta") %>%
#   filter(system == "TB") %>%
#   filter(season == "summer")
# 
# maxZS <- max(abs(centered$zscore), na.rm = TRUE)
# minZS <- min(centered$zscore, na.rm = TRUE)

SiteXSpeciesFull <- CleanHauls %>%
  pivot_wider(id_cols = Reference:systemZone,
              names_from = Scientificname,
              values_from = N2,
              values_fill = 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  ungroup()

LTxSpeciesMeans <- SiteXSpeciesFull %>%
  subset(select = -c(Reference, systemZone, seasonYear)) %>%
  group_by(system, season) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

LTxSpeciesStdev <- SiteXSpeciesFull %>%
  subset(select = -c(Reference, systemZone, seasonYear)) %>%
  group_by(system, season) %>%
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

YearXSpeciesZ <- SiteXSpeciesFull %>%
  subset(select = -c(Reference, systemZone)) %>%
  group_by(system, season, seasonYear) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  pivot_longer(cols = !c(system:seasonYear),
              names_to = "Scientificname",
              values_to = "avg") %>%
              #expand back out to long form for leftjoins
  ungroup() %>%
  left_join(pivot_longer(data = LTxSpeciesMeans,
                         cols = !c(system:season),
                         names_to = "Scientificname",
                         values_to = "LTmean")) %>%
  left_join(pivot_longer(data = LTxSpeciesStdev,
                         cols = !c(system:season),
                         names_to = "Scientificname",
                         values_to = "stdev")) %>%
  group_by(system, season, seasonYear) %>%
  mutate(zscore = ((avg - LTmean)/stdev)) %>%
  filter(LTmean > 0) %>% #remove taxa entirely absent from each system - if LTmean = 0, taxa never observed
  filter(Scientificname != "No fish") %>%
  #filter(Scientificname == "Leiostomus xanthurus") %>%
  filter(system == "AP") %>%
  filter(season == "winter")

maxZS <- max(abs(YearXSpeciesZ$zscore), na.rm = TRUE)

#### plot heatmaps ####
ggplot(YearXSpeciesZ, aes(seasonYear, Scientificname, fill= zscore)) + 
  #facet_wrap(system~season) +
 # scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
  #scale_fill_gradientn(colours = terrain.colors(10))  +
  #scale_fill_viridis(option="magma") +
  scale_fill_gradientn(colours=c("blue","white", "red"), 
                       na.value = "grey98",
                       limits = c(maxZS*-1, maxZS),
                       ) +
  geom_tile()

ggplot(centered, aes(x=seasonYear, 
                     y=zscore)) + 
  #facet_wrap(system~season) +
  # scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
  #scale_fill_gradientn(colours = terrain.colors(10))  +
  #scale_fill_viridis(option="magma") +
  #scale_fill_gradientn(colours=c("blue","white", "red"), na.value = "grey98",
                       #limits = c(-1, 1)) +
  geom_line() +
  facet_wrap(vars(system))

#### centered logic ####

centeredLogic <- centered %>%
  mutate(logic = if_else(zscore == 0, 0,
                         if_else(zscore > 0, 1, -1)))

summaryLogic <- centeredLogic %>%
  group_by(system, season, seasonYear, Scientificname) %>%
  summarise(avg = mean(logic)) %>%
  mutate(logic = if_else(avg == 0, 0,
                         if_else(avg > 0, 1, -1))) %>%
  filter(Scientificname != "No fish") %>%
  filter(system == "TB") %>%
  filter(season == "summer")

ggplot(summaryLogic, aes(seasonYear, Scientificname, fill= logic)) + 
  #facet_wrap(system~season) +
  # scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
  #scale_fill_gradientn(colours = terrain.colors(10))  +
  #scale_fill_viridis(option="magma") +
  scale_fill_gradientn(colours=c("blue","white", "red"), na.value = "grey98",
                       limits = c(-1, 1)) +
  geom_tile()

#### Dornelas style ####
# Linear regression to population abundances
#   ignoring when species was absent (pre/post)
# Sqrt first
# Then scale for mean 0 and stdev 1
# Fit OLS regression through trans data and calc slope/stat sig

library(tseries)

binaryAbs <- rawAbs %>%
  mutate(logic = if_else(N2 > 0, 1, 0)) %>%
  group_by(system, season, seasonYear, Scientificname) %>%
  summarise(avg = mean(logic))

SiteXSpeciesFull <- CleanHauls %>%
  #mutate(N2 = N2^.5) %>% #square root transform
  mutate(logic = if_else(N2 > 0, 1, 0)) %>%
  #group_by(Reference) %>%
  pivot_wider(id_cols = Reference:systemZone,
              #id_expand = TRUE,
              names_from = Scientificname,
              values_from = logic,
              values_fill = 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  ungroup() %>%
  subset(select = -c(Reference, systemZone)) %>%
  group_by(system, season, seasonYear) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  group_by(system, season, seasonYear) %>%
  mutate(across(everything(), ~ if_else(.x > 0, 1, 0))) %>%
  filter(system == "TB") %>%
  filter(season == "summer")

runs.test(as.factor(SiteXSpeciesFull$`Opsanus beta`),
          alternative = "less")