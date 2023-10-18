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

load('TidyGearCode20.Rdata')
load('SXS_filtered.Rdata')

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
           sep = "_") %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))

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

systemSeason_list <- SXS_filtered %>%
  select(systemSeason) %>%
  distinct()

# run simple linear models and check Q-Q plots for diagnostics
linearRegs <- list()
plots <- list()
for(i in systemSeason_list$systemSeason){
  print(i) #watch progress through list

  df <- data.frame()
  df <- SXR_filtered %>%
    filter(systemSeason == i)
  
   #run model
   out <- lm(meanN ~ as.numeric(as.character(seasonYear)), data = df)
  
   #generate simulated residuals
   linearRegs[[i]] <- simulateResiduals(fittedModel = out, plot = F)
   p <- recordPlot() #cludgy way to convert graphics to grob
   plot.new()
   plotQQunif(linearRegs[[i]])
   grid.echo()
   a <- grid.grab() #save grob in loop 
   plots[[i]] <- grid.arrange(a, top = paste0(i)) #label each grob with systemSeason
   
}

#plot em up
wrap_plots(plots, ncol = 2)

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
  geom_ribbon(
    aes(ymin=q10BVC,
        ymax=q90BVC),
    linetype=2, alpha=0.1, color="black") +
  geom_smooth(method="lm", se = FALSE) +
  geom_errorbar(aes(ymin = meanBVC-seBVC, ymax = meanBVC+seBVC))+
  labs(title = "Annual Bottom Vegetation Coverage Over Time",
       x     = "Year",
       y     = "Mean Annual Bottom Vegetation Coverage (%)",
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
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),  geom = 'text', 
                  aes(label = paste("p-value = ", signif(..p.value.., digits = 4), 
                                    "\n R-squared = ", signif(..r.squared.., digits = 3), sep = "")),
                  label.x = 2004, label.y = 85, size = 4) +
  facet_grid(season ~ system)

# ggsave(plot = BVCPlot,
#        filename = "./Outputs/BVCPlot.png",
#        width = 16,
#        height = 9)

TempPlot <- ggplot(waterBVC_full,
       aes(x    = as.numeric(as.character(seasonYear)), 
           y    = meanTemp, 
           #color = season
       )) + 
  geom_point(#size   = 2, 
    #stroke = 0.1,
    #pch    = 21, 
    #colour = "black"
  ) +
  geom_ribbon(
    aes(ymin=q10Temp,
        ymax=q90Temp),
    linetype=2, alpha=0.1, color="black") +
  geom_smooth(method = "lm", se = FALSE) +
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
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),  geom = 'text', 
                  aes(label = paste("p-value = ", signif(..p.value.., digits = 4), 
                                    "\n R-squared = ", signif(..r.squared.., digits = 3), sep = "")),
                  label.x = 2004, label.y = 25, size = 4) +
  facet_grid(season ~ system)

# ggsave(plot = TempPlot,
#        filename = "./Outputs/TempPlot.png",
#        width = 16,
#        height = 9)

NPlot <- ggplot(SXR_filtered,
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
  geom_smooth(method = "lm", se = FALSE) +
  geom_errorbar(aes(ymin = meanN-seN, ymax = meanN+seN))+
  labs(title = "Annual Richness Over Time",
       x     = "Year",
       y     = "Mean Annual Richness per Haul",
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
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),  geom = 'text', 
                  aes(label = paste("p-value = ", signif(..p.value.., digits = 4), 
                                    "\n R-squared = ", signif(..r.squared.., digits = 3), sep = "")),
                  label.x = 2004, label.y = 12, size = 4) +
  facet_grid(season ~ system)

# ggsave(plot = NPlot,
#        filename = "./Outputs/NPlot.png",
#        width = 16,
#        height = 9)

AbundPlot <- ggplot(SXAb,
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
  geom_smooth(method = "lm", se = FALSE) +
  geom_errorbar(aes(ymin = meanAb-seAb, ymax = meanAb+seAb))+
  labs(title = "Annual Abundance Over Time",
       x     = "Year",
       y     = "Mean Annual Total Abundance per Haul",
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
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = y ~ x),  geom = 'text', 
                  aes(label = paste("p-value = ", signif(..p.value.., digits = 4), 
                                    "\n R-squared = ", signif(..r.squared.., digits = 3), sep = "")),
                  label.x = 2004, label.y = 1500, size = 4) +
  facet_grid(season ~ system)

# ggsave(plot = AbundPlot,
#        filename = "./Outputs/AbundPlot.png",
#        width = 16,
#        height = 9)
