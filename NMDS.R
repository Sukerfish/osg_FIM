#### data input ######
library(tidyverse)
library(vegan)
#library(BiodiversityR)
#library(MASS)
library(colortools)
library(ggrepel)
library(RColorBrewer)
library(egg)
library(patchwork)

load('TidyGearCode20.Rdata')

#get the biological data and associated site chars
CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, BottomVegCover, systemZone, Scientificname, N2))
#filter(season == "summer" | season == "winter")
#filter(season == "winter")
#filter(season == "summer")

#get water temp from hydrology dataset
WaterTemp <- HydroList %>%
  subset(select = c(Reference, Temperature))

#collect everything by Reference and spread to long form
SXS_full <- CleanHauls %>%
  mutate(N2 = N2^.25) %>% #fourth-root transform
  group_by(Reference) %>%
  filter(sum(N2)>0) %>% #remove all References with 0 taxa found
  spread(Scientificname,N2) %>%
  ungroup() %>%
  #subset(select = -c(Reference, season, systemZone)) %>%
  replace(is.na(.), 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  group_by(season, system, seasonYear) %>%
  mutate(seasonYear = as.factor(seasonYear)) %>%
  left_join(WaterTemp) %>%
  filter(!is.na(Temperature)) %>% #only excludes 7 sampling events from above pool
  unite("systemSeason",
        system:season,
        sep = "_") %>%
  ungroup()
#summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))


##### setup for CAP loop ####
systemSeason_list <- SXS_full %>%
  select(systemSeason) %>%
  distinct()

NMDSforAll <- list()
for(i in systemSeason_list$systemSeason){
  print(i) #watch progress through list
  
  df <- SXS_full %>% #filter out systemSeason of interest
    filter(systemSeason %in% i)
  
  df_spe <- df %>% #pull out taxa only
    subset(select = -c(systemSeason, seasonYear, Reference, systemZone, BottomVegCover, Temperature)) %>%
    select(which(!colSums(., na.rm=TRUE) %in% 0)) #select only taxa present in this systemSeason

  df_pa <- df_spe
  df_pa[df_pa > 0] <- 1 #convert to pa
  
  spp <- length(df_spe)
  spx <- nrow(df_spe)
  
  df_pa_filtered <- df_pa %>%
    select_if(colSums(.)>(0.05*spx))
  
  df_filtered <- df %>%
    select(c(Reference, seasonYear, all_of(colnames(df_pa_filtered)))) %>%
    rowwise() %>%
    mutate(N = sum(across(!c(Reference:seasonYear)))) %>%
    ungroup() %>%
    filter(N > 0) %>%
    select(!c(N))
  
  df_spe_filtered <- df %>%
    filter(Reference %in% df_filtered$Reference) %>%
    subset(select = -c(systemSeason, seasonYear, Reference, systemZone, BottomVegCover, Temperature)) %>%
    select(which(!colSums(., na.rm=TRUE) %in% 0)) #select only taxa present in this systemSeason
  
  # df_env <- data.frame(df %>% #pull out environmental variables
  #                        subset(select = c(systemSeason, seasonYear, BottomVegCover, Temperature)) %>%
  #                        mutate(contYear = as.numeric(as.character(seasonYear))))
  
  bf <- (vegdist(df_spe_filtered))^0.5 #Bray-Curtis w/ sqrt to reduce negative eigenvalues
  
  nmds <- metaMDS(bf,
                  distance = "bray",
                  k = 2,
                  maxit = 50, 
                  trymax = 50,
                  wascores = TRUE)
  
  NMDSforAll[[i]] <- nmds
  #CAPsforAll[[i]] <- add.spec.scores(CAPsforAll[[i]], df_spe)
  
  #cbPalette1    <- brewer.pal(8, "Paired")[c(2, 3, 6, 7)]
  #cbPalette2    <- brewer.pal(4, "Dark2")
}

data.scores = as.data.frame(scores(nmds))

#data.scores$Sample = pc$Sample
data.scores$Time = df_filtered$seasonYear
#data.scores$Type = pc$Type

cent <- aggregate(cbind(NMDS1, NMDS2) ~ Time, data = data.scores, FUN = mean)

ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 2, aes( colour = as.numeric(as.character(Time))))+ 
  # scale_x_continuous(limits       = c(-.01,.01),
  #                    #breaks       = c(-6,-3,0,3,6),
  #                    minor_breaks = NULL) +
  # scale_y_continuous(limits       = c(-.004,.004),
  #                    #breaks       = c(-6,-3,0,3,6),
  #                    minor_breaks = NULL) +
  scale_colour_gradient(
    low = "red",
    high = "blue",
  ) +
  labs(title = "NMDS",
       x     = "NMDS1",
       y     = "NMDS2",
       caption = paste0("Stress: ", nmds$stress),
  ) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        #legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank(),
        legend.position ="none",
  )

ggplot(cent, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( colour = as.numeric(as.character(Time))))+ 
  # scale_x_continuous(limits       = c(-.01,.01),
  #                    #breaks       = c(-6,-3,0,3,6),
  #                    minor_breaks = NULL) +
  # scale_y_continuous(limits       = c(-.004,.004),
  #                    #breaks       = c(-6,-3,0,3,6),
  #                    minor_breaks = NULL) +
  scale_colour_gradient(
    low = "red",
    high = "blue",
  ) +
  labs(title = "NMDS",
       x     = "NMDS1",
       y     = "NMDS2",
       caption = paste0("Stress: ", nmds$stress),
  ) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        #legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank(),
        legend.position ="none",
        )


#save(CAPsforAll, file = "CAPsforAll.RData")