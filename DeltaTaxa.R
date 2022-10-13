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

summerHauls <- CleanHauls %>%
  filter(season == "summer")

winterHauls <- CleanHauls %>%
  filter(season == "winter")

###### main ######
# SiteXSpeciesFull <- CleanHauls %>%
#   #mutate(N2 = N2^.25) %>% #fourth-root transform
#   group_by(Reference) %>%
#   spread(Scientificname,N2) %>%
#   ungroup() %>%
#   subset(select = -c(Reference, systemZone)) %>%
#   replace(is.na(.), 0) #replace all NA values with 0s, i.e. counting as true zero
#   
# SXSSummary <- SiteXSpeciesFull %>%
#   group_by(system, season, seasonYear) %>%
#   summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

rawAbs <- CleanHauls %>%
  subset(select = c(Reference, system, season, seasonYear, Scientificname, N2))

means <- CleanHauls %>%
  #mutate(N2 = N2^.25) %>% #fourth-root transform
  group_by(system, season, Scientificname) %>%
  summarise(mean = mean(N2))

stdev <- CleanHauls %>%
  group_by(system, season, Scientificname) %>%
  summarise(stdev = sd(N2))

centered <- rawAbs %>%
  left_join(means) %>%
  left_join(stdev) %>%
  mutate(zscore = ((N2 - mean)/stdev))
  # filter(Scientificname == "Anchoa spp.") %>%
  # filter(seasonYear == "2001") %>%
  # filter(system == "AP") %>%
  # filter(season == "winter") %>%
  # summarise(avg = mean(zscore))

summary <- centered %>%
  group_by(system, season, seasonYear, Scientificname) %>%
  summarise(avg = mean(zscore)) %>%
  filter(system == "TB") %>%
  filter(season == "winter")

ggplot(summary, aes(seasonYear, Scientificname, fill= avg)) + 
  #facet_wrap(system~season) +
 # scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
  #scale_fill_gradientn(colours = terrain.colors(10))  +
  #scale_fill_viridis(option="magma") +
  scale_fill_gradientn(colours=c("blue","white", "red"), na.value = "grey98",
                       limits = c(-3.5, 3.5)) +
  geom_tile()

