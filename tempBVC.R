# water temp and bottom veg cover plots

#### data input ######
library(tidyverse)
library(vegan)
#library(BiodiversityR)
#library(MASS)
#library(colortools)
library(ggrepel)
library(RColorBrewer)
library(egg)
library(patchwork)

load('TidyGearCode20.Rdata')
load('SXS_filtered.Rdata')

waterBVC_full <- SXS_filtered %>%
  select(c(Reference, systemSeason, seasonYear, Temperature, BottomVegCover)) %>%
  group_by(systemSeason, seasonYear) %>%
  mutate(meanBVC = mean(BottomVegCover),
         meanTemp = mean(Temperature),
         sdTemp = sd(Temperature),
         sdBVC = sd(BottomVegCover),
         nTemp = length(Temperature),
         nBVC = length(BottomVegCover)) %>%
  mutate(seTemp = sdTemp/sqrt(nTemp),
         seBVC = sdBVC/sqrt(nBVC)) %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_") %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))

ggplot(waterBVC_full,
       aes(x    = as.numeric(as.character(seasonYear)), 
           y    = meanBVC, 
           #color = season
           )) + 
  geom_point(#size   = 2, 
             #stroke = 0.1,
             #pch    = 21, 
             #colour = "black"
               ) +
  geom_smooth(method=lm) +
  geom_errorbar(aes(ymin = meanBVC-seBVC, ymax = meanBVC+seBVC))+
  labs(title = "Annual Bottom Veg Cover Over Time",
       x     = "Year",
       y     = "Mean Annual Bottom Veg Cover (%)",
       #fill  = NULL
       ) +
  # scale_fill_manual(values = cbPalette1,
  #                   labels = c(unique(as.character(df_env$seasonYear)))) +
  theme_bw() +
  theme(legend.text       = element_text(size=rel(0.8)),
        legend.position   = c(0.1,0.89),
        legend.background = element_blank(),
        legend.key        = element_blank(),
        #panel.grid        = element_blank()
        ) +
  facet_grid(season ~ system)

ggplot(waterBVC_full,
       aes(x    = as.numeric(as.character(seasonYear)), 
           y    = meanTemp, 
           #color = season
       )) + 
  geom_point(#size   = 2, 
    #stroke = 0.1,
    #pch    = 21, 
    #colour = "black"
  ) +
  geom_smooth(method=lm) +
  geom_errorbar(aes(ymin = meanTemp-seTemp, ymax = meanTemp+seTemp))+
  labs(title = "Annual Water Temperature Over Time",
       x     = "Year",
       y     = "Mean Annual Water Temperature (°C)",
       #fill  = NULL
  ) +
  # scale_fill_manual(values = cbPalette1,
  #                   labels = c(unique(as.character(df_env$seasonYear)))) +
  theme_bw() +
  theme(legend.text       = element_text(size=rel(0.8)),
        legend.position   = c(0.1,0.89),
        legend.background = element_blank(),
        legend.key        = element_blank(),
        #panel.grid        = element_blank()
  ) +
  facet_grid(season ~ system)

####### linear models #######
