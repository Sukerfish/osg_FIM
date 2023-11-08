library(tidyverse)
library(patchwork)
library(ggpmisc) #for labeling
library(broom) #for labeling
library(DHARMa)
library(gridExtra)
library(gridGraphics)
library(grid)

load('TidyGearCode20.Rdata')
load('SXS_filtered.Rdata')

#### diagnostics ####

SXR_diagnostic <- SXS_filtered %>%
  select(!c(systemSeason, seasonYear, systemZone, BottomVegCover, Temperature)) %>%
  pivot_longer(cols = !c(Reference),
               names_to = "taxa") %>% #pivot to reference and taxa only
  mutate(nRaw = value^4) %>%
  group_by(Reference) %>%
  summarise(abundRaw = sum(nRaw)) %>% #sum all abundance values 
  mutate(abund = abundRaw^0.25) %>%
  left_join(SXS_filtered_env) %>%
  left_join(SXR_filtered_spp)

# run simple linear models and check Q-Q plots for diagnostics
#richness first
linearRegs <- list()
plots <- list()
for(i in systemSeason_list$systemSeason){
  print(i) #watch progress through list
  
  df <- data.frame()
  df <- SXR_diagnostic %>%
    filter(systemSeason == i)
  
  #run model
  out <- lm(N ~ as.numeric(as.character(seasonYear)), data = df)
  
  #generate simulated residuals
  linearRegs[[i]] <- simulateResiduals(fittedModel = out, plot = F)
  p <- recordPlot() #cludgy way to convert graphics to grob
  plot.new()
  #plotResiduals(linearRegs[[i]])
  #plotQQunif(linearRegs[[i]])
  grid.echo()
  a <- grid.grab() #save grob in loop 
  plots[[i]] <- grid.arrange(a, top = paste0(i)) #label each grob with systemSeason
  
}

#plot em up
QQNPlot <- wrap_plots(plots, ncol = 2)
# ggsave(plot = QQNPlot,
#        filename = "./Outputs/QQNPlot.png",
#        width = 16,
#        height = 9)

#abundance next
linearRegs <- list()
plots <- list()
for(i in systemSeason_list$systemSeason){
  print(i) #watch progress through list
  
  df <- data.frame()
  df <- SXR_diagnostic %>%
    filter(systemSeason == i)
  
  #run model
  out <- lm(abund ~ as.numeric(as.character(seasonYear)), data = df)
  
  #generate simulated residuals
  linearRegs[[i]] <- simulateResiduals(fittedModel = out, plot = F)
  p <- recordPlot() #cludgy way to convert graphics to grob
  plot.new()
  #plotResiduals(linearRegs[[i]])
  plotQQunif(linearRegs[[i]])
  grid.echo()
  a <- grid.grab() #save grob in loop 
  plots[[i]] <- grid.arrange(a, top = paste0(i)) #label each grob with systemSeason
  
}

#plot em up
QQAbPlot <- wrap_plots(plots, ncol = 2)
# ggsave(plot = QQAbPlot,
#        filename = "./Outputs/QQAbPlot.png",
#        width = 16,
#        height = 9)

#then temp
linearRegs <- list()
plots <- list()
for(i in systemSeason_list$systemSeason){
  print(i) #watch progress through list
  
  df <- data.frame()
  df <- SXR_diagnostic %>%
    filter(systemSeason == i)
  
  #run model
  out <- lm(Temperature ~ as.numeric(as.character(seasonYear)), data = df)
  
  #generate simulated residuals
  linearRegs[[i]] <- simulateResiduals(fittedModel = out, plot = F)
  p <- recordPlot() #cludgy way to convert graphics to grob
  plot.new()
  #plot(linearRegs[[i]])
  plotQQunif(linearRegs[[i]])
  grid.echo()
  a <- grid.grab() #save grob in loop 
  plots[[i]] <- grid.arrange(a, top = paste0(i)) #label each grob with systemSeason
  
}

#plot em up
QQTempPlot <- wrap_plots(plots, ncol = 2)

#then temp
linearRegs <- list()
plots <- list()
for(i in systemSeason_list$systemSeason){
  print(i) #watch progress through list
  
  df <- data.frame()
  df <- SXR_diagnostic %>%
    filter(systemSeason == i)
  
  #run model
  out <- lm(BottomVegCover ~ as.numeric(as.character(seasonYear)), data = df)
  
  #generate simulated residuals
  linearRegs[[i]] <- simulateResiduals(fittedModel = out, plot = F)
  p <- recordPlot() #cludgy way to convert graphics to grob
  plot.new()
  #plot(linearRegs[[i]])
  plotQQunif(linearRegs[[i]])
  grid.echo()
  a <- grid.grab() #save grob in loop 
  plots[[i]] <- grid.arrange(a, top = paste0(i)) #label each grob with systemSeason
  
}

#plot em up
QQBVCPlot <- wrap_plots(plots, ncol = 2)
