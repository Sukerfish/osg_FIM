#### data input ######

library(viridis)
library(tidyverse)
library(ggplot2)
library(ggallin)
library(writexl)
library(egg)
library(betapart)
library(vegan)

load('TidyGearCode20.Rdata')
load('SXS_filtered.Rdata') #already 4th root transformed

osg_theme <- readRDS('osg_theme.rds')

#load in base data

bpart <- SXS_filtered %>%
  #filter(systemSeason == "TB_summer") %>%
  group_by(systemSeason, seasonYear) %>%
  mutate(seasonYear = as.numeric(as.character(seasonYear))) %>%
  arrange(seasonYear, .by_group = TRUE) %>%
  select(!c(Reference, systemZone, BottomVegCover, Temperature)) %>%
  #mutate_at(vars(-group_cols()), list(~.^.25))
  summarise_all(list(mean)) %>%
  ungroup()

years <- bpart %>%
  select(systemSeason, seasonYear) %>%
  group_by(systemSeason) %>%
  arrange(seasonYear, .by_group = TRUE) %>%
  mutate(yr = as.numeric(row_number()))

           
betaTimeBal <- list()
betaTimeGra <- list()
betaTime <- list()
for (i in bpart$systemSeason){
  #init
  sf <- data.frame()
  sf <- bpart %>%
    filter(systemSeason == i) %>%
    arrange(seasonYear) %>%
    select(!c(systemSeason, seasonYear)) %>%
    select(where(~sum(.) != 0))
  
  out <- beta.pair.abund(sf, index.family = "bray")
  
  betaTimeBal[[i]] <- as.data.frame(as.matrix(out$beta.bray.bal)) %>%
    select(1) %>%
    rename(Balanced = 1) %>%
    mutate(yr = as.numeric(rownames(.)))
  betaTimeGra[[i]] <- as.data.frame(as.matrix(out$beta.bray.gra)) %>%
    select(1) %>%
    rename(Gradient = 1) %>%
    mutate(yr = as.numeric(rownames(.)))
  betaTime[[i]] <- as.data.frame(as.matrix(out$beta.bray)) %>%
    select(1) %>%
    rename(Overall = 1) %>%
    mutate(yr = as.numeric(rownames(.)))
}

#out as dataframe

#assemble components first
betaTimeGraDF <- bind_rows(betaTimeGra, .id = "systemSeason") %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_",
           remove = FALSE) %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH"))) %>%
  left_join(years)

betaTimeBalDF <- bind_rows(betaTimeBal, .id = "systemSeason") %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_",
           remove = FALSE) %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH"))) %>%
  left_join(years)

#assemble overall and attach components
betaTimeDF <- bind_rows(betaTime, .id = "systemSeason") %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_",
           remove = FALSE) %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH"))) %>%
  left_join(years) %>%
  left_join(betaTimeGraDF) %>%
  left_join(betaTimeBalDF) %>%
  select(!c(yr)) %>%
  pivot_longer(cols = c("Overall", "Balanced", "Gradient")) %>% #prep for large plot
  rename(Index = name) %>%
  mutate(season = str_to_title(season)) %>%
  filter(value > 0) %>% #filter index value to clip off initialization year at 0
  mutate(contYear = as.numeric(as.character(seasonYear))) %>%
  ungroup() %>%
  mutate(smallYear = contYear - min(contYear))

system_name <- c(
  AP = "Apalachicola Bay",
  CK = "Cedar Key",
  TB = "Tampa Bay",
  CH = "Charlotte Harbor"
)

betaTimePlot <- ggplot(aes(x = smallYear,
                            y = value,
                            color = Index),
                        data = betaTimeDF) +
  geom_line(
    linewidth = 0.7) +
  facet_grid(season~system) +
  # scale_x_continuous(breaks= seq(1998,2020,5)) +
  # scale_color_cmocean(discrete = TRUE,
  #                     name = "solar") +
  # scale_color_viridis_d(option = "viridis") +
    # geom_smooth(method = "lm", 
    #             formula = y ~ x,
    #             se = FALSE) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Abundance-based dissimilarity over time",
       x = "Year",
       y = "Betadiversity index",
       color = "Component") +
  theme_bw() +
  theme(legend.position="bottom")

betaWinter <- betaTimeDF %>%
  filter(season == "Winter") %>%
  filter(Index == "Overall")
nestWinter <- betaTimeDF %>%
  filter(season == "Winter") %>%
  filter(Index == "Gradient")
turnWinter <- betaTimeDF %>%
  filter(season == "Winter") %>%
  filter(Index == "Balanced")

betaSummer <- betaTimeDF %>%
  filter(season == "Summer") %>%
  filter(Index == "Overall")
nestSummer <- betaTimeDF %>%
  filter(season == "Summer") %>%
  filter(Index == "Gradient")
turnSummer <- betaTimeDF %>%
  filter(season == "Summer") %>%
  filter(Index == "Balanced")

betaWinterLm <- aov(value ~ smallYear * system, data = betaWinter)
nestWinterLm <- aov(value ~ smallYear * system, data = nestWinter)
turnWinterLm <- aov(value ~ smallYear * system, data = turnWinter)
betaSummerLm <- aov(value ~ smallYear * system, data = betaSummer)
nestSummerLm <- aov(value ~ smallYear * system, data = nestSummer)
turnSummerLm <- aov(value ~ smallYear * system, data = turnSummer)

# ggsave(plot = betaTimePlot,
#        filename = "./Outputs/betaTime.png",
#        width = 16,
#        height = 9)

######