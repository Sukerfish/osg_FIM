# new things

library(readxl)
library(tidyverse)
library(ggplot2)
library(viridis)
library(ggallin)
library(egg)

osg_theme <- readRDS('osg_theme.rds')

systems <- c("AP", "CK", "TB", "CH")

sigSlopes <- list()
for (i in unique(systems)){
  sigSlopes[[i]] <- read_excel("Outputs/SigSlopes.xlsx",
                          sheet = i)
}

sigSlopesdf <- bind_rows(sigSlopes)

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
  ungroup()

#collapse full site x species matrix to long term means
LTxSpeciesMeans <- SiteXSpeciesFull %>%
  subset(select = -c(Reference, systemZone, seasonYear)) %>%
  group_by(system, season) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

#make matching standard deviation matrix to long term means
LTxSpeciesStdev <- SiteXSpeciesFull %>%
  subset(select = -c(Reference, systemZone, seasonYear)) %>%
  group_by(system, season) %>%
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

#Z score conversion process
YearXSpeciesZ <- SiteXSpeciesFull %>%
  subset(select = -c(Reference, systemZone)) %>%
  group_by(system, season, seasonYear) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>% #collapse to annual means 
  ungroup() %>%
  pivot_longer(cols = !c(system:seasonYear),
               names_to = "Scientificname",
               values_to = "avg") %>%
  #expand back out to long form for leftjoins with LT mean and LT stdev
  ungroup() %>%
  left_join(pivot_longer(data = LTxSpeciesMeans, #LT mean column added
                         cols = !c(system:season),
                         names_to = "Scientificname",
                         values_to = "LTmean")) %>%
  left_join(pivot_longer(data = LTxSpeciesStdev, #LT stdev column added
                         cols = !c(system:season),
                         names_to = "Scientificname",
                         values_to = "stdev")) %>%
  group_by(system, season, seasonYear) %>%
  mutate(zscore = ((avg - LTmean)/stdev)) %>% #calculate zscores using annual means
  ungroup() %>%
  filter(LTmean > 0) %>% #remove taxa entirely absent from each system - if LTmean = 0, taxa never observed
  filter(Scientificname != "No fish") 
# #filter(Scientificname == "Leiostomus xanthurus") %>%
# filter(system == "AP") %>%
# filter(season == "winter")

YearXSigSpeciesZ <- sigSlopesdf %>%
  left_join(YearXSpeciesZ) %>%
  filter(sign(slope) == -1)

system_name <- c(
  AP = "Apalachicola Bay",
  CK = "Cedar Key",
  TB = "Tampa Bay",
  CH = "Charlotte Harbor"
)

sigSlopePlot <- ggplot(data = YearXSigSpeciesZ) +
  geom_line(aes(x = seasonYear,
                y = zscore,
                color = Scientificname),
            #linewidth = 0.7
            ) +
  facet_grid(season~system,
             labeller = labeller(system = system_name)) +
  scale_x_continuous(breaks = seq(1998,2020,2)) +
  osg_theme +
  theme(legend.position = "none")
  
  # scale_color_cmocean(discrete = TRUE,
  #                     name = "solar") +
  # scale_color_viridis_d(option = "viridis") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Abundance-based dissimilarity over time",
       x = "Year",
       y = "Betadiversity index",
       color = "Component") +
  osg_theme
