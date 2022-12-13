# Heatmap of zscored abundance data
#
#
#### data input ######

library(viridis)
library(tidyverse)
library(ggplot2)
#library(writexl)

load('TidyGearCode20.Rdata')

#### Z scored abundance ####

#grab basic haul data and counts
CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, systemZone, Scientificname, N2))

#full site by species matrix
SiteXSpeciesFull <- CleanHauls %>%
  pivot_wider(id_cols = Reference:systemZone,
              names_from = Scientificname,
              values_from = N2,
              values_fill = 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  ungroup() %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))

totalAbundance <- SiteXSpeciesFull %>%
  rowwise() %>%
  mutate(Ntotal = sum(across(!c(Reference:systemZone)))) %>%
  select(Reference:systemZone, Ntotal)

XtotalAbundance <- totalAbundance %>%
  mutate(xform = Ntotal^.5)

ggplot(totalAbundance, 
       aes(Ntotal))+
  geom_boxplot() +
  facet_grid(season~system) +
  coord_cartesian(xlim = c(0, 500)) +
  #xlab("Population Change") +
  #ylab("Number of Taxa") +
  theme(axis.text=element_text(size = 12)) +
  theme(axis.title=element_text(size = 16)) +
  theme(strip.text = element_text(size = 16)) +
  #ggtitle("GLM w/ no Xform") +
  theme(title=element_text(size = 20))

ggplot(totalAbundance, aes(y=Ntotal,
                    x = factor(seasonYear))) +
  geom_boxplot() +
  scale_x_discrete(breaks=c(2000,2005,2010,2015,2020))+
  coord_cartesian(ylim = c(0, 500)) +
  theme(axis.text=element_text(size = 12)) +
  theme(axis.title=element_text(size = 16)) +
  theme(strip.text = element_text(size = 16)) +
  ggtitle("Abundance per haul over time") +
  theme(title=element_text(size = 20)) +
  facet_grid(season~system)

ggplot(XtotalAbundance, 
       aes(xform))+
  geom_boxplot() +
  facet_grid(season~system) +
  coord_cartesian(xlim = c(0, 30)) +
  #xlab("Population Change") +
  #ylab("Number of Taxa") +
  theme(axis.text=element_text(size = 12)) +
  theme(axis.title=element_text(size = 16)) +
  theme(strip.text = element_text(size = 16)) +
  #ggtitle("GLM w/ no Xform") +
  theme(title=element_text(size = 20))

ggplot(XtotalAbundance, aes(y=xform,
                           x = factor(seasonYear))) +
  geom_boxplot() +
  scale_x_discrete(breaks=c(2000,2005,2010,2015,2020))+
  coord_cartesian(ylim = c(0, 30)) +
  theme(axis.text=element_text(size = 12)) +
  theme(axis.title=element_text(size = 16)) +
  theme(strip.text = element_text(size = 16)) +
  ggtitle("Abundance per haul over time") +
  theme(title=element_text(size = 20)) +
  facet_grid(season~system)

# test<- SiteXSpeciesFull %>%
#   filter(Reference == "APM2000120701") %>%
#   select(!c(Reference:systemZone)) %>%
#   mutate(tot = sum(across(everything())))
