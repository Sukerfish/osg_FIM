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
library(ggpmisc) #for labeling
library(broom) #for labeling
library(DHARMa)
library(gridExtra)
library(gridGraphics)
library(grid)
library(broom)
library(pixiedust)

load('TidyGearCode20.Rdata')
#load('SXS_filtered.Rdata')
load("SXS_filtered_fars.Rdata")

#### data wrangling ####
waterBVC_full <- SXS_filtered %>%
  select(c(Reference, systemSeason, seasonYear, Temperature, BottomVegCover)) %>%
  group_by(systemSeason, seasonYear) %>%
  summarise(meanBVC = mean(BottomVegCover),
         meanTemp = mean(Temperature),
         sdTemp = sd(Temperature),
         sdBVC = sd(BottomVegCover),
         nTemp = length(Temperature),
         nBVC = length(BottomVegCover),
         q10Temp = quantile(Temperature, 0.1),
         q10BVC = quantile(BottomVegCover, 0.1),
         q90Temp = quantile(Temperature, 0.9),
         q90BVC = quantile(BottomVegCover, 0.9)) %>%
  mutate(seTemp = sdTemp/sqrt(nTemp),
         seBVC = sdBVC/sqrt(nBVC)) %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_",
           remove = FALSE) %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))

sdCheck <- waterBVC_full %>%
  ungroup() %>%
  group_by(systemSeason) %>%
  mutate(ltmTemp = mean(meanTemp),
         ltmSDT = sd(meanTemp),
         sdCheck = abs((ltmTemp-meanTemp)/ltmSDT)) %>%
  select(seasonYear, ltmTemp, ltmSDT, sdCheck) %>%
  arrange(sdCheck)

SXAb <- SXS_filtered %>%
  select(!c(systemSeason, seasonYear, systemZone, BottomVegCover, Temperature)) %>%
  pivot_longer(cols = !c(Reference),
               names_to = "taxa") %>% #pivot to reference and taxa only
  mutate(nRaw = value^4) %>%
  group_by(Reference) %>%
  summarise(abundRaw = sum(nRaw)) %>% #sum all abundance values 
  mutate(abund = abundRaw^0.25) %>%
  left_join(SXS_filtered) %>%
  ungroup() %>%
  group_by(systemSeason, seasonYear) %>%
  summarise(meanAb = mean(abundRaw),
            sdAb = sd(abundRaw),
            nAb = length(abundRaw),
            q10Ab = quantile(abundRaw, 0.1),
            q90Ab = quantile(abundRaw, 0.9)) %>%
  mutate(seAb = sdAb/sqrt(nAb)) %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_",
           remove = FALSE) %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))

SXR_filtered_spp <- SXS_filtered %>%
  select(!c(systemSeason, seasonYear, systemZone, BottomVegCover, Temperature)) %>%
  pivot_longer(cols = !c(Reference),
               names_to = "taxa") %>% #pivot to reference and taxa only
  mutate(PA = ifelse(value > 0, 1, 0)) %>% #convert to presence/absence
  select(!c(taxa, value)) %>% #drop taxa and raw counts
  group_by(Reference) %>% #group only by ref in prep
  mutate(N = sum(PA)) %>% #sum all PA values for sample richness
  select(!c(PA)) %>% #drop PA column
  distinct() #keep distinct pairs of ref/richness

#build month dataframe
monthly <- HydroList %>%
  select(Reference, month)

#rejoin with environmental data
SXR_filtered <- SXR_filtered_spp %>%
  left_join(SXS_filtered_env) %>%
  left_join(monthly) %>%
  group_by(systemSeason) %>%
  add_count(name = "n_hauls") %>% #count number of hauls per systemSeason
  ungroup() %>%
  group_by(systemSeason, seasonYear) %>%
  summarise(meanBVC = mean(BottomVegCover),
            meanTemp = mean(Temperature),
            meanN = mean(N),
            sdTemp = sd(Temperature),
            sdBVC = sd(BottomVegCover),
            sdN = sd(N),
            nTemp = length(Temperature),
            nBVC = length(BottomVegCover),
            nN = length(N),
            q10Temp = quantile(Temperature, 0.1),
            q10BVC = quantile(BottomVegCover, 0.1),
            q10N = quantile(N, 0.1),
            q90Temp = quantile(Temperature, 0.9),
            q90BVC = quantile(BottomVegCover, 0.9),
            q90N = quantile(N, 0.9)) %>%
  mutate(seTemp = sdTemp/sqrt(nTemp),
         seBVC = sdBVC/sqrt(nBVC),
         seN = sdN/sqrt(nN)) %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_",
           remove = FALSE) %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))

#rejoin with environmental data
SXF_filtered <- SXR_filtered_spp %>%
  left_join(SXS_filtered_env) %>%
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
         se_mon = sd_mon/sqrt(n_temp_mon)) %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_") 

abundanceModelDF <- SXS_filtered %>%
  dplyr::select(!c(systemSeason, seasonYear, systemZone, BottomVegCover, Temperature)) %>%
  pivot_longer(cols = !c(Reference),
               names_to = "taxa") %>% #pivot to reference and taxa only
  mutate(nRaw = value^4) %>%
  group_by(Reference) %>%
  summarise(abundRaw = sum(nRaw),
            abundAdd = sum(value)) %>% #sum all abundance values 
  mutate(abund = abundRaw^0.25)

#factor date setup with explicit standard interval of months
years <- c(1998:2020)
months <- c(1:12)
dateSetup <- data.frame(contYear = rep(years, each = 12),
                        month = rep(months, length(years)))
dateSetup <- dateSetup %>%
  unite(yearMonth, c(month, contYear), sep = "/", remove = FALSE) %>%
  mutate(dateMonth = as.Date(paste(month, "01", contYear, sep="/"), format="%m/%d/%Y"))

totAbModelDF <- abundanceModelDF %>%
  left_join(SXF_filtered) %>% #merge in enviro data
  # separate(systemSeason,
  #          c("system","season"),
  #          sep = "_") %>%
  mutate(contYear = as.numeric(as.character(seasonYear))) %>%
  unite(yearMonth, c(month, contYear), sep = "/", remove = FALSE) %>%
  unite(systemSeason, c(system, season), sep = "_", remove = FALSE) %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH"))) %>%
  mutate(seasonYear = as.factor(seasonYear)) %>%
  #declare factor levels in prep for ar1 covariance
  mutate(yearMonth = factor(yearMonth, levels = dateSetup$yearMonth)) %>%
  mutate(cYear = contYear - mean(contYear))

systemSeason_list <- SXS_filtered %>%
  select(systemSeason) %>%
  distinct()

#### main plots ####
BVCPlot <- ggplot(data = waterBVC_full,
                  aes(x    = as.numeric(as.character(seasonYear)), 
                      y    = meanBVC, 
                      #color = season
                  )) + 
  geom_point(#size   = 2, 
    #stroke = 0.1,
    #pch    = 21, 
    #colour = "black"
  ) +
  # geom_ribbon(
  #   aes(ymin=q10BVC,
  #       ymax=q90BVC),
  #   linetype=2, alpha=0.1, color="black") +
  #geom_smooth(method="lm", se = FALSE) +
  geom_errorbar(aes(ymin = meanBVC-seBVC, ymax = meanBVC+seBVC))+
  labs(#title = "Annual Bottom Vegetation Coverage Over Time",
       x     = "Year",
       y     = "Mean annual bottom vegetation coverage (%)",
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
  # stat_fit_glance(method = 'lm',
  #                 method.args = list(formula = y ~ x),  geom = 'text', 
  #                 aes(label = paste("p-value = ", signif(after_stat(p.value), digits = 3), 
  #                                   "\n R-squared = ", signif(after_stat(r.squared), digits = 2), sep = "")),
  #                 label.x = 2005, label.y = 20, size = 3) +
  facet_grid(season ~ system)

# ggsave(plot = BVCPlot,
#        filename = "./Outputs/BVCPlot.png",
#        width = 16,
#        height = 9)

TempPlot <- ggplot(waterBVC_full,
                   aes(x    = as.numeric(as.character(seasonYear)), 
                       y    = meanTemp, 
                       #color = system
                   )) + 
  geom_point(#size   = 2, 
    #stroke = 0.1,
    #pch    = 21, 
    #colour = "black"
  ) +
  # geom_ribbon(
  #   aes(ymin=q10Temp,
  #       ymax=q90Temp),
  #   linetype=2, alpha=0.1, color="black") +
  #geom_smooth(method = "lm", se = FALSE) +
  geom_errorbar(aes(ymin = meanTemp-seTemp, ymax = meanTemp+seTemp))+
  labs(#title = "Annual Water Temperature Over Time",
       x     = "Year",
       y     = "Mean annual water temperature (°C)",
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
  # stat_fit_glance(method = 'lm',
  #                 method.args = list(formula = y ~ x),  geom = 'text', 
  #                 aes(label = paste("p-value = ", signif(after_stat(p.value), digits = 3), 
  #                                   "\n R-squared = ", signif(after_stat(r.squared), digits = 2), sep = "")),
  #                 label.x = 2005, label.y = 25, size = 3) +
  facet_grid(season ~ system)

# ggsave(plot = TempPlot,
#        filename = "./Outputs/TempPlot.png",
#        width = 16,
#        height = 9)

NPlot <- ggplot(totAbModelDF,
                    aes(x    = as.numeric(as.character(seasonYear)), 
                        y    = N, 
                        #color = season
                    )) + 
  geom_point(#size   = 2, 
    #stroke = 0.1,
    #pch    = 21, 
    #colour = "black"
  ) +
  # geom_ribbon(
  #   aes(ymin=q10Temp,
  #       ymax=q90Temp),
  #   linetype=2, alpha=0.1, color="black") +
  #geom_smooth(method = "lm", se = FALSE) +
  #geom_errorbar(aes(ymin = meanN-seN, ymax = meanN+seN))+
  labs(#title = "Annual Richness Over Time",
       x     = "Year",
       y     = "Mean annual richness per haul",
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
  # stat_fit_glance(method = 'lm',
  #                 method.args = list(formula = y ~ x),  geom = 'text', 
  #                 aes(label = paste("p-value = ", signif(after_stat(p.value), digits = 3), 
  #                                   "\n R-squared = ", signif(after_stat(r.squared), digits = 2), sep = "")),
  #                 label.x = 2005, label.y = 4, size = 3) +
  facet_grid(season ~ system)

NPlot_ann <- ggplot(SXR_filtered,
                aes(x    = as.numeric(as.character(seasonYear)), 
                    y    = meanN, 
                    #color = season
                )) + 
  geom_point(#size   = 2, 
    #stroke = 0.1,
    #pch    = 21, 
    #colour = "black"
  ) +
  # geom_ribbon(
  #   aes(ymin=q10Temp,
  #       ymax=q90Temp),
  #   linetype=2, alpha=0.1, color="black") +
  # geom_smooth(method = "lm", se = FALSE) +
  geom_errorbar(aes(ymin = meanN-seN, ymax = meanN+seN))+
  labs(#title = "Annual Richness Over Time",
       x     = "Year",
       y     = "Mean annual richness per haul",
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
  # stat_fit_glance(method = 'lm',
  #                 method.args = list(formula = y ~ x),  geom = 'text', 
  #                 aes(label = paste("p-value = ", signif(after_stat(p.value), digits = 3), 
  #                                   "\n R-squared = ", signif(after_stat(r.squared), digits = 2), sep = "")),
  #                 label.x = 2005, label.y = 4, size = 3) +
  facet_grid(season ~ system)

# ggsave(plot = NPlot,
#        filename = "./Outputs/NPlot.png",
#        width = 16,
#        height = 9)

AbundPlot_ann <- ggplot(SXAb,
                    aes(x    = as.numeric(as.character(seasonYear)), 
                        y    = meanAb, 
                        #color = season
                    )) + 
  geom_point(#size   = 2, 
    #stroke = 0.1,
    #pch    = 21, 
    #colour = "black"
  ) +
  # geom_ribbon(
  #   aes(ymin=q10Temp,
  #       ymax=q90Temp),
  #   linetype=2, alpha=0.1, color="black") +
  # geom_smooth(method = "lm", se = FALSE) +
  geom_errorbar(aes(ymin = meanAb-seAb, ymax = meanAb+seAb))+
  labs(#title = "Annual Abundance Over Time",
       x     = "Year",
       y     = "Mean annual total abundance per haul",
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
  # stat_fit_glance(method = 'lm',
  #                 method.args = list(formula = y ~ x),  geom = 'text', 
  #                 aes(label = paste("p-value = ", signif(after_stat(p.value), digits = 3), 
  #                                   "\n R-squared = ", signif(after_stat(r.squared), digits = 2), sep = "")),
  #                 label.x = 2005, label.y = 1500, size = 3) +
  facet_grid(season ~ system)

# ggsave(plot = AbundPlot,
#        filename = "./Outputs/AbundPlot.png",
#        width = 16,
#        height = 9)

#### linear models for temp/bvc ####

#pull out temp and drop winter 2010
tempSubset <- waterBVC_full %>%
  mutate(contYear = as.numeric(as.character(seasonYear))) %>%
  filter(!(contYear == 2010 & season == "winter")) %>%
  ungroup() %>%
  mutate(smallYear = contYear - min(contYear)) %>%
  mutate(region = ifelse(system %in% c("AP", "CK"), "North", "South"))

BVCSubset <- waterBVC_full %>%
  mutate(contYear = as.numeric(as.character(seasonYear))) %>%
  #filter(!(contYear == 2010 & season == "winter")) %>%
  ungroup() %>%
  mutate(smallYear = contYear - min(contYear)) %>%
  mutate(region = ifelse(system %in% c("AP", "CK"), "North", "South"))

#ancova for initial test of temp over time by system and season
TempYearAnn <- aov(meanTemp ~ smallYear * system * season, data = tempSubset)

#BVCYearAnn <- aov(meanBVC ~ smallYear * system * season, data = tempSubset)

#because interaction in system:season ... see which ones.... almost all summer:summer insig except CH/CK summer... all winter sig diff than summer
TempYearAnnCheck <- TukeyHSD(TempYearAnn, "system:season")

#now separate models for each season
winterTemp <- tempSubset %>%
  filter(season == "winter")
winterBVC <- BVCSubset %>%
  filter(season == "winter")

summerTemp <- tempSubset %>%
  filter(season == "summer")
summerBVC <- BVCSubset %>%
  filter(season == "summer")

#separated winter
winterLm <- aov(meanTemp ~ smallYear * system, data = winterTemp)
winterBVCLm <- aov(meanBVC ~ smallYear * system, data = winterBVC)
summerLm <- aov(meanTemp ~ smallYear * system, data = summerTemp)
summerBVCLm <- aov(meanBVC ~ smallYear * system, data = summerBVC)

#combine temps and go by month now for sinusoidal plot
combinedTemp <- totAbModelDF %>%
  group_by(yearMonth, system) %>%
  summarise(monthTemp = mean(Temperature)) %>%
  mutate(tinyYear = round(decimal_date(as.Date(parse_date_time(yearMonth, "mY"))), 3)) %>%
  ungroup() %>%
  mutate(tinyYear = tinyYear - min(tinyYear))

combinedLm <- lm(monthTemp ~ tinyYear * system + cos(2*pi*tinyYear) + sin(2*pi*tinyYear), data = combinedTemp)

#### full temperature model setup ####
load('GearCode20Refresh.Rdata')

#establish years of interest
YearFilter <- ZoneFilter %>%
  select(system, StartYear, EndYear) %>%
  unique()

#establish seasons of interest
SOI <- c(
  "spring",
  "summer",
  "fall",
  "winter"
)

#expand sampling date and attach years of interest
bigHydroList <- TidyHydro %>%
  mutate(month = month(Sampling_Date), year = year(Sampling_Date)) %>%
  left_join(YearFilter) %>%
  #label specific months as seasons
  mutate(season = if_else(month %in% c(4:5), "spring",
                          if_else(month %in% c(6:9), "summer",
                                  if_else(month %in% c(10:11), "fall",
                                          "winter")))) %>%
  #tag everything with year associated with sampling season, with December
  #as the preceding year's data
  mutate(seasonYear = if_else(month %in% c(12), year + 1, year)) %>%
  #remove the incomplete final seasonYear of winter data
  filter(season != "winter" | seasonYear != EndYear + 1) %>%
  #filter by years of interest (all begin before, so first year of interest is
  #not truncated)
  filter(seasonYear >= StartYear) %>%
  #remove year of interest logic
  subset(select = -c(StartYear, EndYear)) %>%
  #filter by seasons of interest
  filter(season %in% SOI) 

allTemps <- bigHydroList %>%
  mutate(contYear = as.numeric(as.character(seasonYear))) %>%
  unite(yearMonth, c(month, contYear), sep = "/", remove = FALSE) %>%
  unite(systemSeason, c(system, season), sep = "_", remove = FALSE) %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH"))) %>%
  mutate(seasonYear = as.factor(seasonYear)) %>%
  #declare factor levels
  mutate(yearMonth = factor(yearMonth, levels = dateSetup$yearMonth))%>%
  group_by(yearMonth, system) %>%
  summarise(monthTemp = mean(Temperature)) %>%
  mutate(tinyYear = round(decimal_date(as.Date(parse_date_time(yearMonth, "mY"))), 3)) %>%
  ungroup()
  
(allTempsPlot <- ggplot(allTemps,
                aes(x     = tinyYear, 
                    y     = monthTemp, 
                    group = system,
                    color = system
                )) + 
    geom_point() +
    geom_smooth(method = "lm", 
                formula = y ~ x + cos(2*pi*x) + sin(2*pi*x),
                se = FALSE) +
    #geom_errorbar(aes(ymin = meanTemp-seTemp, ymax = meanTemp+seTemp))+
    labs(#title = "Annual Water Temperature Over Time",
      x     = "Year",
      y     = "Mean monthly water temperature (°C)",
      #fill  = NULL
    ) +
    #scale_color_viridis_d() +
    theme_bw() +
  facet_grid(~system) +
    theme(
      legend.position = "bottom"
    )
)


#### old methodology ####
# 
# tempOutputs <- list()
# cleanTemps <- list()
# for (i in systemSeason_list$systemSeason){
#   tempDF <- data.frame()
#   tempDF <- filter(waterBVC_full, systemSeason == i) %>%
#     mutate(contYear = as.numeric(as.character(seasonYear))) %>%
#     filter(!(contYear == 2010 & season == "winter")) %>%
#     mutate(cYear = contYear - mean(contYear))
#   
#   tempOut <- lm(meanTemp ~ cYear, tempDF)
#   tempOutputs[[i]] <- tempOut
#   cleanTemps[[i]] <- broom::tidy(tempOut)
#   
# }
# 
# tempTest <- bind_rows(cleanTemps, .id = "systemSeason")
# redOut <- dust(tempTest) %>%
#   sprinkle(cols = c("estimate", "std.error", "statistic"), round = 3) %>%
#   sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>%
#   sprinkle_colnames(term = "Term", p.value = "P-value")
# 
# tempEstimates <- as.data.frame(redOut) %>%
#   select(c(systemSeason, Term, "estimate")) %>%
#   pivot_wider(names_from = Term, values_from = "estimate")
# 
# tempPvalues <- as.data.frame(redOut) %>%
#   select(c(systemSeason, Term, "P-value")) %>%
#   pivot_wider(names_from = Term, values_from = "P-value")
# 
# BVCOutputs <- list()
# cleanBVCs <- list()
# for (i in systemSeason_list$systemSeason){
#   bvcDF <- data.frame()
#   bvcDF <- filter(waterBVC_full, systemSeason == i) %>%
#     mutate(contYear = as.numeric(as.character(seasonYear))) %>%
#     #filter(!(contYear == 2010 & season == "winter")) %>%
#     mutate(cYear = contYear - mean(contYear))
#   
#   bvcOut <- lm(meanBVC ~ cYear, bvcDF)
#   BVCOutputs[[i]] <- bvcOut
#   cleanBVCs[[i]] <- broom::tidy(bvcOut)
#   
# }
# 
# bvcTest <- bind_rows(cleanBVCs, .id = "systemSeason")
# greenOut <- dust(bvcTest) %>%
#   sprinkle(cols = c("estimate", "std.error", "statistic"), round = 3) %>%
#   sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>%
#   sprinkle_colnames(term = "Term", p.value = "P-value")
# 
# bvcEstimates <- as.data.frame(greenOut) %>%
#   select(c(systemSeason, Term, "estimate")) %>%
#   pivot_wider(names_from = Term, values_from = "estimate")
# 
# bvcPvalues <- as.data.frame(greenOut) %>%
#   select(c(systemSeason, Term, "P-value")) %>%
#   pivot_wider(names_from = Term, values_from = "P-value")