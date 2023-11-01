##### Import and Filter Data #####

library(tidyverse)
library(patchwork)
library(lubridate)

load('TidyGearCode20.Rdata')
load('SXS_filtered.Rdata')

osg_theme <- readRDS('osg_theme.rds')

#load in base data

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
         sd_mon = sd(Temperature)) %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_") 

filteredTemp <- SXR_filtered %>%
  dplyr::select(!c(Reference, N, systemZone, BottomVegCover, Temperature, n_hauls)) %>%
  distinct() %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))%>%
  mutate(season = str_to_title(season))

##### temperature plots #######
lowerLimits <- ggplot(filteredTemp,
       aes(x    = as.numeric(as.character(seasonYear)), 
           y    = lower, 
           #color = season
       )) + 
  geom_point(#size   = 2, 
    #stroke = 0.1,
    #pch    = 21, 
    #colour = "black"
  ) +
  geom_smooth(method=lm) +
  geom_errorbar(aes(ymin = avg_temp-se_temp, ymax = avg_temp+se_temp))+
  labs(title = "Coldest 10% Water Temperature Over Time",
       x     = "Year",
       y     = "Annual 10% Coldest Water Temperature  Threshold (°C)",
       #fill  = NULL
  ) +
  osg_theme +
  facet_grid(season ~ system)

upperLimits <- ggplot(filteredTemp,
                      aes(x    = as.numeric(as.character(seasonYear)), 
                          y    = upper, 
                          #color = season
                      )) + 
  geom_point(#size   = 2, 
    #stroke = 0.1,
    #pch    = 21, 
    #colour = "black"
  ) +
  geom_smooth(method=lm) +
  geom_errorbar(aes(ymin = avg_temp-se_temp, ymax = avg_temp+se_temp))+
  labs(title = "Warmest 10% Water Temperature Over Time",
       x     = "Year",
       y     = "Annual 10% Warmest Water Temperature  Threshold (°C)",
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

longTemp <- filteredTemp %>%
  group_by(system, seasonYear) %>%
  pivot_longer(cols = c("lower", "upper")) %>%
  mutate(name = factor(name, levels = c("upper", "lower")))


capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

system_name <- c(
  AP = "Apalachicola Bay",
  CK = "Cedar Key",
  TB = "Tampa Bay",
  CH = "Charlotte Harbor"
)
#full plot
total <- ggplot(longTemp,
                      aes(x    = as.numeric(as.character(seasonYear)), 
                          y    = value, 
                          #shape = season,
                          color = name
                      )) + 
  geom_point(size = 2
  ) +
  geom_hline(aes(yintercept = avg_ltm),
             linetype = "dashed") +
  geom_smooth(method=lm) +
  geom_errorbar(aes(ymin = value-se_temp, ymax = value+se_temp))+
  scale_x_continuous(breaks= seq(1998,2020,4)) +
  scale_color_brewer(palette = "Set1",
                     labels = c("upper" = "Upper 90th", "lower" = "Lower 10th")) +
  labs(title = "Water Temperature Extremes Over Time",
       x     = "Year",
       y     = "Annual Water Temperature (°C)",
       color = "Percentile",
       #shape = "Season"
  ) +
  osg_theme +
  facet_grid(season~system,
             labeller = labeller(.default = capitalize,
                                 system = system_name))

# ggsave(plot = total,
#        filename = "./Outputs/tempTime.png",
#        width = 16,
#        height = 9)


tempX <- upperLimits + lowerLimits

sdTemp <- ggplot(data = filteredTemp,
  # pivot_longer(filteredTemp, 
  #              cols = !c(system, season, seasonYear), 
  #              names_to = "metric") %>%
  #   filter(metric %in% c("lower_sd", "upper_sd", "sd_ann")),
                 aes(x    = as.numeric(as.character(seasonYear)), 
                     y    = sd_ann, 
                     #color = metric
                 )) + 
  geom_point(#size   = 2, 
    #stroke = 0.1,
    #pch    = 21, 
    #colour = "black"
  ) +
  geom_smooth(method=lm) +
  #geom_errorbar(aes(ymin = avg_temp-se_temp, ymax = avg_temp+se_temp))+
  labs(title = "Standard Deviation Water Temperature Over Time",
       x     = "Year",
       y     = "SD Water Temperature (°C)",
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


#### model data setup #####

modelDFM <- SXR_filtered %>%
  # separate(systemSeason,
  #          c("system","season"),
  #          sep = "_") %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH"))) %>%
  mutate(contYear = as.numeric(as.character(seasonYear))) %>%
  mutate(seasonYear = as.factor(seasonYear)) %>%
  # mutate(n_hauls = log(n_hauls)) %>%
  # filter(season == "summer") %>%
  # filter(system == "TB") %>%
  filter(season == "winter")

modelDF <- as.data.frame(modelDFM)

# abundanceModelDF <- SXS_filtered %>%
#   select(!c(systemSeason, seasonYear, BottomVegCover, systemZone, Temperature)) %>%
#   #group_by(Reference) %>%
#   rowwise() %>% #sum abundance across all taxa
#   mutate(abund = sum(across(-Reference))) %>%
#   select(c(Reference, abund))

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
  unite(yearMonth, c(month, contYear), sep = "/", remove = FALSE)

totAbModelDF <- abundanceModelDF %>%
  left_join(SXR_filtered) %>% #merge in enviro data
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

abundDist <- ggplot(data = totAbModelDF) +
  geom_histogram(aes(x = abundAdd),
                 binwidth = 1) +
  facet_grid(system~season)

# modelDFAb_summer <- as.data.frame(totAbModelDF %>%
#                                     filter(season == "summer"))
#   # unite(yearMonth, c(contYear, month), remove = FALSE) %>%
#   # mutate(yearMonth = as.factor(yearMonth))
# modelDFAb_winter <- as.data.frame(totAbModelDF %>%
#                                     filter(season == "winter"))
  # unite(yearMonth, c(contYear, month), remove = FALSE) %>%
  # mutate(yearMonth = as.factor(yearMonth))

####plots
# ggplot(modelDF, aes(x=N)) +
#   geom_histogram(binwidth=1) +
#   theme(axis.text=element_text(size = 12)) +
#   theme(axis.title=element_text(size = 16)) +
#   theme(strip.text = element_text(size = 16)) +
#   ggtitle("Richness per haul") +
#   theme(title=element_text(size = 20)) +
#   facet_grid(season~system)
# 
# ggplot(modelDF, aes(y=N,
#                     x = factor(seasonYear))) +
#   geom_boxplot() +
#   scale_x_discrete(breaks=c(1995,2000,2005,2010,2015,2020))+
#   theme(axis.text=element_text(size = 12)) +
#   theme(axis.title=element_text(size = 16)) +
#   theme(strip.text = element_text(size = 16)) +
#   ggtitle("Richness per haul over time") +
#   theme(title=element_text(size = 20)) +
#   facet_grid(season~system)
# 
# ggplot(modelDFAb, aes(x=abund)) +
#   geom_histogram(binwidth=1) +
#   theme(axis.text=element_text(size = 12)) +
#   theme(axis.title=element_text(size = 16)) +
#   theme(strip.text = element_text(size = 16)) +
#   ggtitle("Abundance per haul") +
#   theme(title=element_text(size = 20)) +
#   facet_grid(season~system)
# 
# ggplot(modelDFAb, aes(y=abund,
#                       x = factor(contYear))) +
#   geom_boxplot() +
#   scale_x_discrete(breaks=c(2000,2005,2010,2015,2020))+
#   theme(axis.text=element_text(size = 12)) +
#   theme(axis.title=element_text(size = 16)) +
#   theme(strip.text = element_text(size = 16)) +
#   ggtitle("Abundance per haul over time") +
#   theme(title=element_text(size = 20)) +
#   facet_grid(season~system)


##########################################################
# GLMM:                                                  #
# Response: richness or total abundance                  #
# Fixed effects: estuary (four levels) and year          #
# Random effect: within-estuary zone                     #
# Offset term for variation in number of seine hauls     #
# Autoregressive term to account for serial correlation  #
##########################################################
#### modeling ####

# library(nlme)
# library(lme4)
# library(MASS) # needs MASS (version 7.3-58)
library(glmmTMB)
library(DHARMa)
library(gridExtra)
library(gridGraphics)
library(grid)
library(AICcmodavg)
library(performance)
library(jtools)
#library(fitdistrplus)

library(broom)
library(broom.mixed)
library(pixiedust)

sysList <- unique(totAbModelDF$systemSeason)

#### richness ####
linearRegs <- list()
plots <- list()
richnessOut <- list()
outputs <- list()
sysGLMMS_richness <- list()
for (i in sysList){
  runner <- data.frame()
  runner <- filter(totAbModelDF, systemSeason == i)
 
  glmmOut <- glmmTMB(N ~ year_Z +
                           temp_Z +
                           bvc_Z +
                       #(1|contYear) +
                       #(1|systemZone),
                       offset(log(n_hauls)) +
                           ar1(yearMonth + 0|systemZone),
                         data = runner,
                     #dispformula = ~0,
                     #ziformula = ~1,
                         family = gaussian)
  sysGLMMS_richness[[i]] <- glmmOut
  outputs[[i]] <- broom.mixed::tidy(glmmOut, effects = "fixed")
  # richnessOut[[i]] <- dust(glmmOut, effects = "fixed", caption = i) %>%
  #   sprinkle(cols = c("estimate", "std.error", "statistic"), round = 3) %>%
  #   sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>%
  #   sprinkle_colnames(term = "Term", p.value = "P-value")

  #generate simulated residuals
  # linearRegs[[i]] <- simulateResiduals(fittedModel = glmmOut, plot = F)
  # p <- recordPlot() #cludgy way to convert graphics to grob
  # plot.new()
  # #plotResiduals(linearRegs[[i]])
  # plotQQunif(linearRegs[[i]])
  # grid.echo()
  # a <- grid.grab() #save grob in loop
  # plots[[i]] <- grid.arrange(a, top = paste0(i)) #label each grob with systemSeason
}

richOut <- bind_rows(outputs, .id = "systemSeason")
pinkOut <- dust(richOut) %>%
  sprinkle(cols = c("estimate", "std.error", "statistic"), round = 3) %>%
  sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>%
  sprinkle_colnames(term = "Term", p.value = "P-value")

funOut <- as.data.frame(pinkOut) %>%
  select(c(systemSeason, Term, "estimate")) %>%
  pivot_wider(names_from = Term, values_from = "estimate")

# write.csv(funOut, file = "./Outputs/richnessOut.csv",
#           row.names = FALSE)

#wrap_plots(plots, ncol = 2)

#check_model(sysGLMMS_richness$CH_winte)
#effect_plot(glmmOut, pred = contYear, interval = TRUE, partial.residuals = TRUE)

#summary(sysGLMMS_richness$CH_winte)

##### abundance #####

linearRegs <- list()
plots <- list()
abundOut <- list()
outputs <- list()
sysGLMMS_abund <- list()
for (i in sysList){
  funner <- data.frame()
  funner <- filter(totAbModelDF, systemSeason == i)
  
  glmmOut <- glmmTMB(abund ~ year_Z +
                       temp_Z +
                       bvc_Z +
                       offset(log(n_hauls)) +
                       ar1(yearMonth + 0|systemZone),
                     data = funner,
                     family = gaussian)
  sysGLMMS_abund[[i]] <- glmmOut
  outputs[[i]] <- broom.mixed::tidy(glmmOut, effects = "fixed")
  
  # abundOut[[i]] <- dust(glmmOut, effects = "fixed", caption = i) %>%
  #   sprinkle(cols = c("estimate", "std.error", "statistic"), round = 3) %>% 
  #   sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>% 
  #   sprinkle_colnames(term = "Term", p.value = "P-value")
  
  # #generate simulated residuals
  # linearRegs[[i]] <- simulateResiduals(fittedModel = glmmOut, plot = F)
  # p <- recordPlot() #cludgy way to convert graphics to grob
  # plot.new()
  # plot(linearRegs[[i]])
  # #plotQQunif(linearRegs[[i]])
  # grid.echo()
  # a <- grid.grab() #save grob in loop 
  # plots[[i]] <- grid.arrange(a, top = paste0(i)) #label each grob with systemSeason
}


abundOut <- bind_rows(outputs, .id = "systemSeason")
blueOut <- dust(abundOut) %>%
  sprinkle(cols = c("estimate", "std.error", "statistic"), round = 3) %>%
  sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>%
  sprinkle_colnames(term = "Term", p.value = "P-value")

bigOut <- as.data.frame(blueOut) %>%
  select(c(systemSeason, Term, "estimate")) %>%
  pivot_wider(names_from = Term, values_from = "estimate")

# write.csv(bigOut, file = "./Outputs/abundOut.csv",
#           row.names = FALSE)

#wrap_plots(plots, ncol = 2)
#summary(sysGLMMS_abund$CH_winter)
check_model(sysGLMMS_abund$CH_winter)

dust(sysGLMMS_abund$CH_winter, effects = "fixed") %>% 
  sprinkle(cols = c("estimate", "std.error", "statistic"), round = 3) %>% 
  sprinkle(cols = "p.value", fn = quote(pvalString(value)))

#summary(sysGLMMS_abund$CH_winter)



######### temperature AIC #######
library(buildmer)

temp <- totAbModelDF %>%
  filter(season == "winter")

tempList <- unique(temp$systemSeason)
#varList <- c("upper_Z", "lower_Z", "sd_t_Z")

linearRegs <- list()
plots <- list()
tempAbundOut <- list()
outputs <- list()
tempGLMMS_abund <- list()
for (i in tempList){
  print(i)
  funner <- data.frame()
  funner <- filter(totAbModelDF, systemSeason == i)
  
  glmmOut <- buildglmmTMB(abund ~ year_Z +
                       temp_Z +
                       bvc_Z +
                         upper_Z +
                         lower_Z +
                         offset(log(n_hauls)) +
                       ar1(yearMonth + 0|systemZone),
                     data = funner,
                     family = gaussian(),
                     buildmerControl(crit = "AIC"))
  tempGLMMS_abund[[i]] <- glmmOut
  tempAbundOut[[i]] <- broom.mixed::tidy(glmmOut@model, effects = "fixed")
  
  # abundOut[[i]] <- dust(glmmOut, effects = "fixed", caption = i) %>%
  #   sprinkle(cols = c("estimate", "std.error", "statistic"), round = 3) %>% 
  #   sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>% 
  #   sprinkle_colnames(term = "Term", p.value = "P-value")
  
  # #generate simulated residuals
  # linearRegs[[i]] <- simulateResiduals(fittedModel = glmmOut, plot = F)
  # p <- recordPlot() #cludgy way to convert graphics to grob
  # plot.new()
  # plot(linearRegs[[i]])
  # #plotQQunif(linearRegs[[i]])
  # grid.echo()
  # a <- grid.grab() #save grob in loop 
  # plots[[i]] <- grid.arrange(a, top = paste0(i)) #label each grob with systemSeason
}
summary(tempGLMMS_abund$CH_winter)

tempOut <- bind_rows(tempAbundOut, .id = "systemSeason")
tempAOut <- dust(tempOut) %>%
  sprinkle(cols = c("estimate", "std.error", "statistic"), round = 3) %>%
  sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>%
  sprinkle_colnames(term = "Term", p.value = "P-value")

bigAOut <- as.data.frame(tempAOut) %>%
  select(c(systemSeason, Term, "estimate")) %>%
  pivot_wider(names_from = Term, values_from = "estimate")

# write.csv(bigOut, file = "./Outputs/abundOut.csv",
#           row.names = FALSE)


linearRegs <- list()
plots <- list()
tempNOut <- list()
outputs <- list()
tempGLMMS_abund <- list()
for (i in tempList){
  print(i)
  funner <- data.frame()
  funner <- filter(totAbModelDF, systemSeason == i)
  
  glmmOut <- buildglmmTMB(N ~ year_Z +
                            temp_Z +
                            bvc_Z +
                            upper_Z +
                            lower_Z +
                            offset(log(n_hauls)) +
                            ar1(yearMonth + 0|systemZone),
                          data = funner,
                          family = gaussian(),
                          buildmerControl(crit = "AIC"))
  tempGLMMS_abund[[i]] <- glmmOut
  tempNOut[[i]] <- broom.mixed::tidy(glmmOut@model, effects = "fixed")
  
  # abundOut[[i]] <- dust(glmmOut, effects = "fixed", caption = i) %>%
  #   sprinkle(cols = c("estimate", "std.error", "statistic"), round = 3) %>% 
  #   sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>% 
  #   sprinkle_colnames(term = "Term", p.value = "P-value")
  
  # #generate simulated residuals
  # linearRegs[[i]] <- simulateResiduals(fittedModel = glmmOut, plot = F)
  # p <- recordPlot() #cludgy way to convert graphics to grob
  # plot.new()
  # plot(linearRegs[[i]])
  # #plotQQunif(linearRegs[[i]])
  # grid.echo()
  # a <- grid.grab() #save grob in loop 
  # plots[[i]] <- grid.arrange(a, top = paste0(i)) #label each grob with systemSeason
}

tempNTOut <- bind_rows(tempNOut, .id = "systemSeason")
tempBOut <- dust(tempNTOut) %>%
  sprinkle(cols = c("estimate", "std.error", "statistic"), round = 3) %>%
  sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>%
  sprinkle_colnames(term = "Term", p.value = "P-value")

bigBOut <- as.data.frame(tempBOut) %>%
  select(c(systemSeason, Term, "estimate")) %>%
  pivot_wider(names_from = Term, values_from = "estimate")







#### separated by season ####
 #####abundance 


linearRegs <- list()
plots <- list()
sysGLMMS_abund_summer <- list()
for (i in sysList){
  runner <- data.frame()
  runner <- filter(modelDFAb_summer, system == i)
  
  glmmOut <- glmmTMB(abund ~ year_Z +
                       offset(log(n_hauls)) +
                       temp_Z +
                       bvc_Z +
                       ar1(yearMonth + 0|systemZone),
                     data = runner,
                     family = gaussian)
  sysGLMMS_abund_summer[[i]] <- glmmOut
  
  #generate simulated residuals
  linearRegs[[i]] <- simulateResiduals(fittedModel = glmmOut, plot = F)
  p <- recordPlot() #cludgy way to convert graphics to grob
  plot.new()
  plot(linearRegs[[i]])
  #plotQQunif(linearRegs[[i]])
  grid.echo()
  a <- grid.grab() #save grob in loop 
  plots[[i]] <- grid.arrange(a, top = paste0(i)) #label each grob with systemSeason
}
wrap_plots(plots, ncol = 2)

check_model(sysGLMMS_abund_summer[[i]])

linearRegs <- list()
plots <- list()
sysGLMMS_abund_winter <- list()
for (i in sysList){
  runner <- data.frame()
  runner <- filter(modelDFAb_winter, system == i)
  
  glmmOut <- glmmTMB(abund ~ year_Z +
                       offset(log(n_hauls)) +
                       temp_Z +
                       bvc_Z +
                       ar1(yearMonth + 0|systemZone),
                     data = runner,
                     family = gaussian)
  sysGLMMS_abund_winter[[i]] <- glmmOut
  
  #generate simulated residuals
  linearRegs[[i]] <- simulateResiduals(fittedModel = glmmOut, plot = F)
  p <- recordPlot() #cludgy way to convert graphics to grob
  plot.new()
  plot(linearRegs[[i]])
  #plotQQunif(linearRegs[[i]])
  grid.echo()
  a <- grid.grab() #save grob in loop 
  plots[[i]] <- grid.arrange(a, top = paste0(i)) #label each grob with systemSeason
}
wrap_plots(plots, ncol = 2)

check_model(sysGLMMS_abund_winter[[i]])



#### all systems ####
#richness first
#summer
tmbR_summer <- glmmTMB(N ~ system +
                     cYear +
                     offset(log(n_hauls)) +
                       temp_C +
                       bvc_C +
                     ar1(seasonYear + 0|systemZone),
                   data = modelDFAb_summer,
                   family = gaussian)
simulateResiduals(fittedModel = tmbR_summer, plot = T)
check_model(tmbR_summer)
summary(tmbR_summer)

#winter
tmbR_winter <- glmmTMB(N ~ system +
                         contYear +
                         offset(log(n_hauls)) +
                         temp_Z +
                         bvc_Z +
                         ar1(seasonYear + 0|systemZone),
                       data = modelDFAb_winter,
                       family = gaussian)
simulateResiduals(fittedModel = tmbR_winter, plot = T)
summary(tmbR_winter)
#plot(tmbR_winter)

#then abundance
#summer
tmbA_summer <- glmmTMB(abund ~ system +
                         contYear +
                         offset(log(n_hauls)) +
                         temp_Z +
                         bvc_Z +
                         ar1(seasonYear + 0|systemZone),
                       data = modelDFAb_summer,
                       family = gaussian)
simulateResiduals(fittedModel = tmbA_summer, plot = T)
summary(tmbA_summer)
#plot(tmbA_summer)

#winter
tmbA_winter <- glmmTMB(abund ~ system +
                         contYear +
                         offset(log(n_hauls)) +
                         temp_Z +
                         bvc_Z +
                         ar1(seasonYear + 0|systemZone),
                       data = modelDFAb_winter,
                       family = gaussian)
simulateResiduals(fittedModel = tmbA_winter, plot = T)
summary(tmbA_winter)
#plot(tmbA_winter)

basicSummaries <- list()
basicSummaries$tmbA_winter <- summary(tmbA_winter)
basicSummaries$tmbA_summer <- summary(tmbA_summer)
basicSummaries$tmbR_winter <- summary(tmbR_winter)
basicSummaries$tmbR_summer <- summary(tmbR_summer)

# temperature component testing
#summer
tmbR_tempX_winter <- glmmTMB(N ~ system +
                         contYear +
                         offset(log(n_hauls)) +
                         #Temperature +
                           upper +
                           lower +
                         (1|BottomVegCover) +
                         ar1(seasonYear + 0|systemZone),
                       data = modelDFAb_winter,
                       family = gaussian)
summary(tmbR_tempX_winter)
#winter
tmbA_tempX_winter <- glmmTMB(abund ~ system +
                         contYear +
                         offset(log(n_hauls)) +
                         #Temperature +
                           upper +
                           lower +
                         (1|BottomVegCover) +
                         ar1(seasonYear + 0|systemZone),
                       data = modelDFAb_winter,
                       family = gaussian)
summary(tmbA_tempX_winter)


#richness
#alternative cor structure ... corAR1 converges
glmmPQLcorAR1 <- glmmPQL(N ~ system + contYear + Temperature + offset(log(n_hauls)),
                   random = ~ 1|as.factor(systemZone),
                   family = poisson,
                   #correlation = corAR1(form = ~1|as.factor(systemZone)),
                   data = modelDFAb)
summary(glmmPQLcorAR1)
plot(glmmPQLcorAR1)



#stepwiseAIC
test1 <- glm(abund ~ system + contYear + upper + lower + sd_ann,
                         #random = ~ 1|as.factor(systemZone),
                         offset = log(n_hauls),
                         family = gaussian,
                         #correlation = corAR1(form = ~1|as.factor(systemZone)),
                         data = modelDFAb)
test <- stepAIC(test1, trace = TRUE)
summary(test)

##### OLD MODELS ####
# #total abundance
# #alternative cor structure ... corAR1 converges
# glmmPQLcorAR1_totAb <- glmmPQL(abund ~ system + contYear + offset(log(n_hauls)),
#                          random = ~ 1|as.factor(systemZone),
#                          family = poisson,
#                          correlation = corAR1(form = ~1|as.factor(systemZone)),
#                          data = modelDFAb)
# summary(glmmPQLcorAR1_totAb)
# plot(glmmPQLcorAR1_totAb)
# 
# 
# 
# 
# #converges without correlation
# glmmPQL <- glmmPQL(n ~ system + seasonYear + offset(log(n_hauls)),
#                    random = ~ 1|systemZone,
#                    family = poisson,
#                    data = modelDF)
# summary(glmmPQL)
# plot(glmmPQL)
# 
# testlmer <- glmer(N ~ system + seasonYear + offset(log(n_hauls))
#                    + (1|as.factor(modelDF$systemZone)),
#                    family = poisson,
#                    #correlation = corAR1(form = ~1|systemZone),
#                    data = modelDF)
# summary(testlmer)
# 
# 
# 
# glmmTMB
# library(glmmTMB)
# rmod_tmb <- glmmTMB(N~system+
#                       seasonYear+
#                       offset(log(n_hauls))+
#                       #ar1(factor(seasonYear))+
#                       (1|systemZone)+
#                       (1|Temperature)+
#                       (1|BottomVegCover),
#                     zi=~0,
#                     family=poisson,
#                     data=modelDF)
# summary(rmod_tmb)
# VarCorr(rmod_tmb)
# 
# lme4 
# #converges if seasonYear is random
# library(lme4)
# library(lmerTest)
# 
# glme_rich <- glmer(N~system+
#                        seasonYear+
#                        offset(log(n_hauls))+
#                        (1|BottomVegCover)+
#                        (1|systemZone)+
#                      (1|Temperature)+
#                        (1|contYear),
#                      family="poisson",
#                      data=modelDF,
#                      control=glmerControl(optimizer="bobyqa",
#                                           check.conv.grad=.makeCC("warning",0.05)))
# 
# summary(glme_rich)
# #ranef(glme_rich)
# 
# # WHEN contYear is included:
# # boundary (singular) fit: see help('isSingular')
# # 
# # > ranef(gmod_lme4_L)
# # $BottomVegCover
# # (Intercept)
# # 0   -0.286778197
# # 1   -0.554656884
# # 2   -0.485282983
# # 3   -0.316112389
# # 4   -0.153897614
# # 5   -0.420456838
# # 6    0.052852022
# # 7   -0.095967185
# # 8   -0.041994018
# # 9    0.182975981
# # 10  -0.277217398
# # 12  -0.333955075
# # 13   0.199986612
# # 15  -0.141748694
# # 17  -0.225760348
# # 20  -0.094252424
# # 24  -0.338278760
# # 25  -0.019703572
# # 28   0.063680466
# # 30   0.038523397
# # 35   0.044137427
# # 38  -0.004546357
# # 40   0.056430271
# # 45   0.142677940
# # 50   0.139696794
# # 55   0.333678760
# # 60   0.201555754
# # 65   0.118314700
# # 69   0.248448168
# # 70   0.217038225
# # 73   0.039055369
# # 75   0.174176325
# # 78   0.236469920
# # 79   0.037584506
# # 80   0.214267889
# # 85   0.291731016
# # 90   0.229240836
# # 92   0.002031692
# # 94   0.050437252
# # 95   0.198828684
# # 96   0.071597185
# # 97  -0.213756871
# # 98   0.117693689
# # 99   0.381708174
# # 100  0.089781643
# # 101 -0.005640213
# # 
# # $contYear
# # (Intercept)
# # 1998           0
# # 1999           0
# # 2000           0
# # 2001           0
# # 2002           0
# # 2003           0
# # 2004           0
# # 2005           0
# # 2006           0
# # 2007           0
# # 2008           0
# # 2009           0
# # 2010           0
# # 2011           0
# # 2012           0
# # 2013           0
# # 2014           0
# # 2015           0
# # 2016           0
# # 2017           0
# # 2018           0
# # 2019           0
# # 2020           0
# # 
# # $systemZone
# # (Intercept)
# # AP_A  0.125793430
# # AP_B -0.123653260
# # CH_A  0.056074224
# # CH_B  0.011136499
# # CH_C  0.047730151
# # CH_D -0.110828474
# # CK_B -0.021993325
# # CK_C  0.024001857
# # TB_A  0.105474887
# # TB_B -0.046134171
# # TB_C -0.019111862
# # TB_D -0.043769500
# # TB_E  0.008875772
# # 
# # with conditional variances for “BottomVegCover” “contYear” “systemZone” 
# summary(glme_rich_summer)
# plot(glme_rich_summer)
# plot(glme_rich_summer,systemZone~resid(.))
# 
# 
# glme_rich_summer <- glmer(n~system+
#                        seasonYear+
#                        offset(log(n_hauls))+
#                        (1|BottomVegCover)+
#                        (1|systemZone),
#                      family="poisson",
#                      data=modelDF_summer,
#                      control=glmerControl(optimizer="bobyqa",
#                                           check.conv.grad=.makeCC("warning",0.05)))
# ranef(glme_rich_summer)
# summary(glme_rich_summer)
# # Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# # Family: poisson  ( log )
# # Formula: n ~ system + seasonYear + offset(log(n_hauls)) + (1 | BottomVegCover) +  
# #   (1 | systemZone)
# # Data: modelDF_summer
# # Control: glmerControl(optimizer = "bobyqa", check.conv.grad = .makeCC("warning",      0.05))
# # 
# # AIC      BIC   logLik deviance df.resid 
# # 53437.9  53636.9 -26691.0  53381.9     8979 
# # 
# # Scaled residuals: 
# #   Min      1Q  Median      3Q     Max 
# # -3.2559 -1.0306 -0.1251  0.8969  6.7312 
# # 
# # Random effects:
# #   Groups         Name        Variance Std.Dev.
# # BottomVegCover (Intercept) 0.067345 0.25951 
# # systemZone     (Intercept) 0.005563 0.07458 
# # Number of obs: 9007, groups:  BottomVegCover, 46; systemZone, 13
# # 
# # Fixed effects:
# #   Estimate Std. Error z value Pr(>|z|)    
# # (Intercept)    -1.714599   0.073560 -23.309  < 2e-16 ***
# #   systemCK       -0.072798   0.075572  -0.963    0.335    
# # systemTB       -0.660035   0.063317 -10.424  < 2e-16 ***
# #   systemCH       -0.666871   0.065484 -10.184  < 2e-16 ***
# #   seasonYear1999  0.010660   0.033726   0.316    0.752    
# # seasonYear2000 -0.002186   0.033806  -0.065    0.948    
# # seasonYear2001 -0.168063   0.029977  -5.606 2.07e-08 ***
# #   seasonYear2002 -0.178254   0.030073  -5.927 3.08e-09 ***
# #   seasonYear2003 -0.138177   0.029717  -4.650 3.32e-06 ***
# #   seasonYear2004 -0.215828   0.029101  -7.416 1.20e-13 ***
# #   seasonYear2005 -0.424016   0.029437 -14.404  < 2e-16 ***
# #   seasonYear2006 -0.378166   0.029125 -12.984  < 2e-16 ***
# #   seasonYear2007 -0.392300   0.029232 -13.420  < 2e-16 ***
# #   seasonYear2008 -0.338052   0.028908 -11.694  < 2e-16 ***
# #   seasonYear2009 -0.366478   0.029043 -12.618  < 2e-16 ***
# #   seasonYear2010 -0.371903   0.029026 -12.813  < 2e-16 ***
# #   seasonYear2011 -0.373310   0.028999 -12.873  < 2e-16 ***
# #   seasonYear2012 -0.438333   0.029357 -14.931  < 2e-16 ***
# #   seasonYear2013 -0.409770   0.029196 -14.035  < 2e-16 ***
# #   seasonYear2014 -0.368040   0.028956 -12.711  < 2e-16 ***
# #   seasonYear2015 -0.396270   0.029029 -13.651  < 2e-16 ***
# #   seasonYear2016 -0.343117   0.028871 -11.884  < 2e-16 ***
# #   seasonYear2017 -0.350129   0.028993 -12.076  < 2e-16 ***
# #   seasonYear2018 -0.379346   0.029058 -13.055  < 2e-16 ***
# #   seasonYear2019 -0.341449   0.028863 -11.830  < 2e-16 ***
# #   seasonYear2020 -0.355158   0.028927 -12.278  < 2e-16 ***
# #   ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# glme_rich_winter <- glmer(n~system+
#                             seasonYear+
#                             offset(log(n_hauls))+
#                             (1|BottomVegCover)+
#                             (1|systemZone),
#                           family="poisson",
#                           data=modelDF_winter,
#                           control=glmerControl(optimizer="bobyqa",
#                                                check.conv.grad=.makeCC("warning",0.05)))
# ranef(glme_rich_winter)
# summary(glme_rich_winter)
# # Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# # Family: poisson  ( log )
# # Formula: n ~ system + seasonYear + offset(log(n_hauls)) + (1 | BottomVegCover) +  
# #   (1 | systemZone)
# # Data: modelDF_winter
# # Control: glmerControl(optimizer = "bobyqa", check.conv.grad = .makeCC("warning",      0.05))
# # 
# # AIC      BIC   logLik deviance df.resid 
# # 47522.8  47721.7 -23733.4  47466.8     8950 
# # 
# # Scaled residuals: 
# #   Min      1Q  Median      3Q     Max 
# # -2.9616 -1.0542 -0.1665  0.8835  6.9494 
# # 
# # Random effects:
# #   Groups         Name        Variance Std.Dev.
# # BottomVegCover (Intercept) 0.04696  0.2167  
# # systemZone     (Intercept) 0.01988  0.1410  
# # Number of obs: 8978, groups:  BottomVegCover, 47; systemZone, 13
# # 
# # Fixed effects:
# #   Estimate Std. Error z value Pr(>|z|)    
# # (Intercept)    -2.524486   0.112436 -22.453  < 2e-16 ***
# #   systemCK       -0.262329   0.142140  -1.846 0.064954 .  
# # systemTB       -0.424898   0.118906  -3.573 0.000352 ***
# #   systemCH       -0.158834   0.122977  -1.292 0.196503    
# # seasonYear1999  0.004921   0.041893   0.117 0.906493    
# # seasonYear2000  0.129259   0.040912   3.159 0.001581 ** 
# #   seasonYear2001  0.067912   0.038215   1.777 0.075554 .  
# # seasonYear2002 -0.155348   0.038880  -3.996 6.45e-05 ***
# #   seasonYear2003  0.060321   0.037483   1.609 0.107549    
# # seasonYear2004 -0.117618   0.036794  -3.197 0.001390 ** 
# #   seasonYear2005 -0.333318   0.037394  -8.914  < 2e-16 ***
# #   seasonYear2006 -0.344907   0.037255  -9.258  < 2e-16 ***
# #   seasonYear2007 -0.229335   0.036449  -6.292 3.14e-10 ***
# #   seasonYear2008 -0.351031   0.037212  -9.433  < 2e-16 ***
# #   seasonYear2009 -0.454898   0.037774 -12.043  < 2e-16 ***
# #   seasonYear2010 -0.259740   0.036674  -7.082 1.42e-12 ***
# #   seasonYear2011 -0.333853   0.037093  -9.000  < 2e-16 ***
# #   seasonYear2012 -0.398587   0.037512 -10.626  < 2e-16 ***
# #   seasonYear2013 -0.277981   0.036759  -7.562 3.96e-14 ***
# #   seasonYear2014 -0.198911   0.036336  -5.474 4.39e-08 ***
# #   seasonYear2015 -0.269293   0.036613  -7.355 1.91e-13 ***
# #   seasonYear2016 -0.224745   0.036466  -6.163 7.13e-10 ***
# #   seasonYear2017 -0.263274   0.036618  -7.190 6.49e-13 ***
# #   seasonYear2018 -0.231000   0.036502  -6.328 2.48e-10 ***
# #   seasonYear2019 -0.283648   0.036705  -7.728 1.10e-14 ***
# #   seasonYear2020 -0.294739   0.036857  -7.997 1.28e-15 ***
# #   ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# interceptOnly <- lmer(n ~ 1 + (1|systemZone),
#                        data = modelDF)
# summary(interceptOnly)
# # Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# # Formula: n ~ 1 + (1 | systemZone)
# #    Data: modelDF
# # 
# # REML criterion at convergence: 47233.5
# # 
# # Scaled residuals: 
# #     Min      1Q  Median      3Q     Max 
# # -2.3774 -0.7433 -0.1415  0.6357  4.8244 
# # 
# # Random effects:
# #  Groups     Name        Variance Std.Dev.
# #  systemZone (Intercept)  2.498   1.581   
# #  Residual               11.045   3.323   
# # Number of obs: 9002, groups:  systemZone, 13
# # 
# # Fixed effects:
# #             Estimate Std. Error      df t value Pr(>|t|)    
# # (Intercept)   5.5275     0.4398 12.0087   12.57 2.86e-08 ***
# # ---
# # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# model1 <- lmer(n ~ 1 + system + seasonYear + (1|systemZone),
#                data = modelDF)
# summary(model1)
# # Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# # Formula: n ~ 1 + system + seasonYear + (1 | systemZone)
# #    Data: modelDF
# # 
# # REML criterion at convergence: 47224.1
# # 
# # Scaled residuals:
# #     Min      1Q  Median      3Q     Max
# # -2.3921 -0.7530 -0.1404  0.6483  4.8387
# # 
# # Random effects:
# #  Groups     Name        Variance Std.Dev.
# #  systemZone (Intercept)  1.045   1.022
# #  Residual               11.045   3.323
# # Number of obs: 9002, groups:  systemZone, 13
# # 
# # Fixed effects:
# #               Estimate Std. Error         df t value Pr(>|t|)
# # (Intercept) -4.350e+00  1.149e+01  8.906e+03  -0.379    0.705
# # systemCK    -1.297e+00  1.029e+00  8.945e+00  -1.261    0.239
# # systemTB     4.995e-01  8.616e-01  8.979e+00   0.580    0.576
# # systemCH     2.395e+00  8.914e-01  8.961e+00   2.687    0.025 *
# # seasonYear   4.551e-03  5.704e-03  8.991e+03   0.798    0.425
# # ---
# # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# # 
# # Correlation of Fixed Effects:
# #            (Intr) systCK systTB systCH
# # systemCK   -0.045
# # systemTB   -0.060  0.597
# # systemCH   -0.056  0.577  0.689
# # seasonYear -0.998  0.000  0.006  0.004
# 
# model2 <- lmer(n ~ 1 + system + seasonYear + (1 + seasonYear|systemZone),
#                data = modelDF)
# summary(model2)
# 
# 
# 
# 
# 
# 
# library(broom.mixed)
# fit_augmented <- augment(glme_rich_summer)
# 
# p1 <-ggplot(fit_augmented, aes(.fitted, .resid))+
#   geom_point()
# 
# p1 + stat_smooth(method="loess") + 
#   geom_hline(yintercept=0, col="red", linetype="dashed") +
#   xlab("Fitted values") + 
#   ylab("Residuals") + 
#   ggtitle("Residual vs Fitted Plot") + 
#   theme_bw()
# 
# qqnorm(fit_augmented[[".resid"]])
# 
# library(aods3)
# gof(glme_rich_summer)
# 
# sims <- simulate(glme_rich_summer,nsim=1000)
# nzeros <- colSums(sims==0)
# par(las=1,bty="l")
# plot(pt <- prop.table(table(nzeros)),
#      ylab="Probability",xlab="Number of zeros")
# (obszero <- sum(modelDF$n==0))
# 
# #gmod_lme4_agq <- update(gmod_lme4_L,nAGQ=10)
# 
# 
# 
# brms 
# library(brms)
# 
# options(mc.cores = parallel::detectCores())
# 
# brm1 <- brm(formula = n ~ system + 
#               seasonYear + 
#               (1|systemZone) + 
#               offset(log(n_hauls)),
#             data = modelDF,
#             family = "poisson")
# 
# (summ_brm1 <- summary(brm1))
# plot(brm1)

##### summaries ####
#define vars for summarizing
varst = c("N", "BottomVegCover")
RichnessSummary <- SXR_filtered %>%
  select(c(system, seasonYear, season, N, BottomVegCover)) %>%
  group_by(system, seasonYear, season) %>%
  summarise(across(any_of(varst),
                   list(
                     mean = ~mean(., na.rm=TRUE),
                     median = ~median(., na.rm=TRUE),
                     q10 = ~quantile(., probs = 0.1, na.rm=TRUE),
                     q90 = ~quantile(., probs = 0.9, na.rm=TRUE)
                   ), .names = "{col}_{fn}"))

#assign factor levels for plot order (North -> South)
RichnessSummary$system = factor(RichnessSummary$system, levels=c("AP","CK","TB","CH"))

#### Plots ####

library(ggplot2)
library(cowplot)

#Richness Plot
p1 <- ggplot(data=RichnessSummary,
             aes(x=as.numeric(as.character(seasonYear)),
                 y=N_mean)) +
  geom_ribbon(
    aes(ymin=N_q10,
        ymax=N_q90,
        fill="10th-90th Percentile"),
    linetype=2, alpha=0.1, color="purple") +
  geom_point(
    aes(y = N_median, 
        color="median")
  ) +
  #theme(legend.position="bottom") +
  theme_bw() + 
  ggtitle("Richness over time") +
  xlab("Year") +
  #scale_x_continuous(breaks = seq(1998, 2020, 2)) +
  ylab("Number of taxa per haul") +
  geom_line(aes(color = "mean")) +
  facet_grid(season ~ system) +
  theme(strip.text.x = element_text(size = 20)) +
  theme(strip.text.y = element_text(size = 20))

#BottomVegCover plot
p2 <- ggplot(data=RichnessSummary,
             aes(x=as.numeric(as.character(seasonYear)), y=BottomVegCover_mean)) +
  geom_ribbon(
    aes(ymin=BottomVegCover_q10,
        ymax=BottomVegCover_q90,
        fill="10th-90th Percentile"),
    linetype=2, alpha=0.1, color="purple") +
  geom_point(
    aes(y = BottomVegCover_median, 
        color="median")
  ) +
  #theme(legend.position="bottom") +
  theme_bw() + 
  ggtitle("SAV coverage over time") +
  xlab("Year") +
  #scale_x_continuous(breaks = seq(1998, 2020, 2)) +
  ylab("Estimated percent SAV coverage") +
  geom_line(aes(color = "mean")) +
  facet_grid(season ~ system) +
  theme(strip.text.x = element_text(size = 20)) +
  theme(strip.text.y = element_text(size = 20))

richGrass <- wrap_plots(p1, p2,
          ncol = 1,
          guides = "collect")

# ggsave(plot = richGrass,
#        filename = "./Outputs/richGrass.png",
#        width = 16,
#        height = 9)