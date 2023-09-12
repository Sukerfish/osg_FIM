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
  mutate(season = str_to_title(season))

betaTimePlot <- ggplot(data = betaTimeDF) +
  geom_line(aes(x = seasonYear,
    y = value,
    color = Index),
    linewidth = 0.7) +
  facet_grid(system~season) +
  scale_x_continuous(breaks= seq(1998,2020,2)) +
  theme_bw() +
  labs(title = "Abundance-based dissimilarity over time",
       x = "Year",
       y = "Betadiversity index",
       color = "Component") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        title = element_text(size = 20),
        legend.text = element_text(size = 12))

# ggsave(plot = betaTimePlot,
#        filename = "./Outputs/betaTime.png",
#        width = 16,
#        height = 9)

######