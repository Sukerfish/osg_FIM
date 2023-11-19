##### NMDS annual abundance data #####

library(tidyverse)
library(vegan)
library(ggrepel)

load('TidyGearCode20.Rdata')
load('SXS_filtered.RData')

monthly <- HydroList %>%
  select(Reference, month)

#### winter first ####
test_winter <- SXS_filtered %>%
  #left_join(monthly) %>%
  filter(systemSeason %in% c("AP_winter", "CK_winter", "CH_winter", "TB_winter")) %>%
  group_by(systemSeason, seasonYear) %>%
  select(!c(Reference, systemZone)) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) #annual summary of all stuff including continuous env vars

#pull out environmental vars and assign factor labels
test_env_winter <- test_winter %>%
  select(systemSeason, BottomVegCover, Temperature) %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_") %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))

#pull out spp and remove taxa not found in this season
test_spp_winter <- test_winter %>%
  ungroup() %>%
  select(!c(systemSeason, seasonYear, BottomVegCover, Temperature)) %>%
  select(which(!colSums(., na.rm=TRUE) %in% 0))

#run nmds
nmds_ann_winter <- metaMDS(test_spp_winter,
                      distance = "bray",
                      k=2,
                      maxit = 50,
                      trymax = 50,
                      wascores = TRUE)

#calculate environmental fectors
env_winter <- envfit(nmds_ann_winter, test_env_winter, permutations = 999, na.rm = TRUE)

#extract NMDS scores
winterScores <- as.data.frame(scores(nmds_ann_winter)$sites)
winterSpp <- as.data.frame(scores(nmds_ann_winter)$species)

#add system column
winterScores$system <- test_env_winter$system

#pull out coords for env vectors
winter_coords_con <- as.data.frame(scores(env_winter, "vectors"))
winter_coords_cat <- as.data.frame(scores(env_winter, "factors")) %>%
  mutate(names = row.names(.)) %>% 
  filter(str_starts(names, "system")) %>%
  mutate(cleaned = str_remove(names, "system")) #clean up

#build plot
winterNMDS <- ggplot(data = winterScores, 
                     aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = winterScores, #sites
             aes(colour = system), 
             size = 3, 
             alpha = .8) + 
  stat_ellipse(aes(colour = system)) +
  geom_point(data = winterSpp, #spp
             size = 1, 
             alpha = 0.5,
             color = "grey60") + 
  scale_colour_viridis_d() +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), #vecs
               data = winter_coords_con, 
               size =1, 
               alpha = 0.5, 
               colour = "grey30",
               arrow = arrow(length = unit(0.03, "npc"))
               ) +
  # geom_point(data = winter_coords_cat, #centroids
  #            aes(x = NMDS1, 
  #                y = NMDS2),
  #            shape = "diamond", 
  #            size = 4, 
  #            alpha = 0.8, 
  #            colour = "grey10") +
  # geom_label_repel(data = winter_coords_cat, #labels
  #           aes(x = NMDS1, y = NMDS2),
  #           label = winter_coords_cat$cleaned, 
  #           colour = "grey10", 
  #           fontface = "bold") +
  # geom_text_repel(data = winterSpp, #labels
  #                  aes(x = NMDS1, y = NMDS2),
  #                  label = row.names(winterSpp), 
  #                  colour = "grey10", 
  #                  fontface = "bold") +
  geom_label_repel(data = winter_coords_con, #labels
            aes(x = NMDS1, y = NMDS2), 
            colour = "grey30", 
            fontface = "bold", 
            label = row.names(winter_coords_con)) + 
  theme_bw() +
  labs(colour = "System",
       title = "NMDS winter",
       caption = paste0("Stress: ", round(nmds_ann_winter$stress, 4)))

winterNMDS

#older plot checkouts
# plot(nmds_ann_winter, type = "n", las = 1, main = "NMDS winter")
# points(nmds_ann_winter, display = "sites")
# points(nmds_ann_winter, display = "species", col = "red", pch = 3)
# ordiellipse(nmds_ann_winter,
#             groups = test_env_winter$systemSeason,
#             kind = "se",
#             conf = 0.95,
#             display = "sites",
#             label=T)
# 
# plot(nmds_ann_winter, main = "NMDS winter")
# ordihull(nmds_ann_winter,groups=test_env_winter$systemSeason,draw="polygon",col="grey90",label=F)
# #plot(SXSRaw_sppfit, p.max = 0.001, col = "black", cex = 0.7)

#### summer next ####
test_summer <- SXS_filtered %>%
  #left_join(monthly) %>%
  filter(systemSeason %in% c("AP_summer", "CK_summer", "CH_summer", "TB_summer")) %>%
  group_by(systemSeason, seasonYear) %>%
  select(!c(Reference, systemZone)) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

test_env_summer <- test_summer %>%
  select(systemSeason, seasonYear, BottomVegCover, Temperature) %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_") %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))

test_spp_summer <- test_summer %>%
  ungroup() %>%
  select(!c(systemSeason, seasonYear, BottomVegCover, Temperature)) %>%
  select(which(!colSums(., na.rm=TRUE) %in% 0))

nmds_ann_summer <- metaMDS(test_spp_summer,
                          distance = "bray",
                          k=2,
                          maxit = 50,
                          trymax = 50,
                          wascores = TRUE)

env_summer <- envfit(nmds_ann_summer, test_env_summer, permutations = 999, na.rm = TRUE)

#extract NMDS scores
summerScores <- as.data.frame(scores(nmds_ann_summer)$sites)
summerSpp <- as.data.frame(scores(nmds_ann_summer)$species)

#add system column
summerScores$system <- test_env_summer$system

#pull out coords for env vectors
summer_coords_con <- as.data.frame(scores(env_summer, "vectors"))
summer_coords_cat <- as.data.frame(scores(env_summer, "factors")) %>%
  mutate(names = row.names(.)) %>% 
  filter(str_starts(names, "system")) %>%
  mutate(cleaned = str_remove(names, "system")) #clean up

#build plot
summerNMDS <- ggplot(data = summerScores, 
                     aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = summerScores, #sites
             aes(colour = system), 
             size = 3, 
             alpha = .8) + 
  stat_ellipse(aes(colour = system)) +
  geom_point(data = summerSpp, #spp
             size = 1, 
             alpha = 0.5,
             color = "grey60") + 
  scale_colour_viridis_d() +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), #vecs
               data = summer_coords_con, 
               size =1, 
               alpha = 0.5, 
               colour = "grey30",
               arrow = arrow(length = unit(0.03, "npc"))) +
  # geom_point(data = summer_coords_cat, #centroids
  #            aes(x = NMDS1, 
  #                y = NMDS2),
  #            shape = "diamond", 
  #            size = 4, 
  #            alpha = 0.8, 
  #            colour = "grey10") +
  # geom_label_repel(data = summer_coords_cat, #labels
  #                  aes(x = NMDS1, y = NMDS2),
  #                  label = summer_coords_cat$cleaned, 
  #                  colour = "grey10", 
  #                  fontface = "bold") +
  # geom_text_repel(data = summerSpp, #labels
  #                  aes(x = NMDS1, y = NMDS2),
  #                  label = row.names(summerSpp), 
  #                  colour = "grey10", 
  #                  fontface = "bold") +
  geom_label_repel(data = summer_coords_con, #labels
                   aes(x = NMDS1, y = NMDS2), 
                   colour = "grey30", 
                   fontface = "bold", 
                   label = row.names(summer_coords_con)) + 
  theme_bw() +
  labs(colour = "System",
       title = "NMDS summer",
       caption = paste0("Stress: ", round(nmds_ann_summer$stress, 4)))

summerNMDS
