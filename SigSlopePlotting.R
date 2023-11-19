# new things

library(readxl)
library(tidyverse)
library(ggplot2)
library(viridis)
library(ggallin)
library(egg)
library(patchwork)

osg_theme <- readRDS('osg_theme.rds')

sigSlopes <- read_excel("./Outputs/sigSlopes_SXS.xlsx")

ecoGroup <- read.csv("Outputs/taxaList_grouped.csv")

#reclass syngnathids under fish
ecoGroup$syngnathid <- 0
ecoGroup$syngnathid[ecoGroup$group == "syngnathid"] <- 1
ecoGroup$group[ecoGroup$group == "syngnathid"] <- "fish"

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
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH"))) %>%
  filter(LTmean > 0) %>% #remove taxa entirely absent from each system - if LTmean = 0, taxa never observed
  filter(Scientificname != "No fish")
# #filter(Scientificname == "Leiostomus xanthurus") %>%
# filter(system == "AP") %>%
# filter(season == "winter")

YearXSigSpeciesZ <- sigSlopes %>%
  left_join(YearXSpeciesZ) %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH"))) %>%
  group_by(system, season) %>%
  left_join(ecoGroup) %>%
  mutate(group = as.factor(group)) %>%
  #arrange(desc(abs(slope)), .by_group = TRUE) %>%
  #ungroup() %>%
  #group_by(Scientificname) %>%
  #slice_head(n = 5) %>%
  #filter(sign(slope) == 1) %>%
  mutate(season = str_to_title(season))

ecoSigNeg <- YearXSigSpeciesZ %>%
  filter(sign(slope) == -1) %>%
#  left_join(ecoGroup) %>%
  group_by(system, season, seasonYear, group) %>%
  summarise(meanZ = mean(zscore))

ecoSigPos <- YearXSigSpeciesZ %>%
  filter(sign(slope) == 1) %>%
  #left_join(ecoGroup) %>%
  group_by(system, season, seasonYear, group) %>%
  summarise(meanZ = mean(zscore))

groupList <- unique(ecoGroup$group)
names(groupList) <- unique(ecoGroup$group)

system_name <- c(
  AP = "Apalachicola Bay",
  CK = "Cedar Key",
  TB = "Tampa Bay",
  CH = "Charlotte Harbor"
)

negsigSlopePlot <- ggplot(data = YearXSigSpeciesZ) +
  geom_line(aes(x = seasonYear,
                y = zscore,
                color = Scientificname),
            #linewidth = 0.7
            ) +
  facet_grid(season~system,
             labeller = labeller(system = system_name)) +
  scale_x_continuous(breaks = seq(1998,2020,4)) +
  osg_theme +
  theme(legend.position = "none")

negecoSlopePlot <- ggplot(data = ecoSigNeg) +
  geom_line(aes(x = seasonYear,
                y = meanZ,
                color = group),
            linewidth = 0.7
  ) +
  facet_grid(season~system,
             labeller = labeller(system = system_name)) +
  scale_x_continuous(breaks = seq(1998,2020,4)) +
  scale_color_manual( #need to force legends to be identical (i.e. show all possibilities) to merge w/ patchwork
    breaks = groupList,
    values = c("red", "blue", "green3", "black", "orange", "purple", "darkblue", "yellow", "gray"), #use scale shape identities
    drop = FALSE
  ) +
  labs(title = "Negative slopes") +
  theme_bw()

posecoSlopePlot <- ggplot(data = ecoSigPos) +
  geom_line(aes(x = seasonYear,
                y = meanZ,
                color = group),
            linewidth = 0.7
  ) +
  facet_grid(season~system,
             labeller = labeller(system = system_name)) +
  scale_x_continuous(breaks = seq(1998,2020,4)) +
  scale_color_manual( #need to force legends to be identical (i.e. show all possibilities) to merge w/ patchwork
    breaks = groupList,
    values = c("red", "blue", "green3", "black", "orange", "purple", "darkblue", "yellow", "gray"), #use scale shape identities
    drop = FALSE
  ) +
  labs(title = "Positive slopes") +
  theme_bw()
  
ecoSlopes <- wrap_plots(negecoSlopePlot, posecoSlopePlot,
                        guides = "collect")

# ggsave(plot = ecoSlopes,
#        filename = "./Outputs/ecoSlopes.png",
#        width = 16,
#        height = 9)

# scale_color_cmocean(discrete = TRUE,
  #                     name = "solar") +
  # scale_color_viridis_d(option = "viridis") +
  # scale_color_brewer(palette = "Set1") +
  # labs(title = "Abundance-based dissimilarity over time",
  #      x = "Year",
  #      y = "Betadiversity index",
  #      color = "Component") +
  # osg_theme
