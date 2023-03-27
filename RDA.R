# PERMANOVAs for each system/season combination


#### data input ######
library(tidyverse)
library(vegan)
library(egg)
library(ggplot2)
library(ggrepel)
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


##### setup for PERMANOVA loop ####
systemSeason_list <- SXS_full %>%
  select(systemSeason) %>%
  distinct()

RDAsforAll <- list()
for(i in systemSeason_list$systemSeason){
  print(i) #watch progress through list
  
  df <- SXS_full %>% #filter out systemSeason of interest
    filter(systemSeason %in% i)
  
  df_spe <- df %>% #pull out taxa only
    subset(select = -c(systemSeason, seasonYear, Reference, systemZone, BottomVegCover, Temperature)) %>%
    select(which(!colSums(., na.rm=TRUE) %in% 0))
  
  df_env <- df %>% #pull out environmental variables
    subset(select = c(systemSeason, seasonYear, BottomVegCover, Temperature)) %>%
    mutate(contYear = as.numeric(as.character(seasonYear)))
  
  bf <- (vegdist(df_spe))^0.5 #Bray-Curtis w/ sqrt to prevent negative eigenvalues
  
  rda = dbrda(bf ~ Temperature + BottomVegCover,
                 data = df_env,
                 #strata = df_env$seasonYear,
                 #add = "lingoes",
                 #parallel = 6,
                 #method="bray", 
                 #permutations = 999,
  )
  env <- envfit(rda, df_env[,3:4])
  scores <- scores(rda, display = "all")
  RDAsforAll[[i]]$rda <- rda
  RDAsforAll[[i]]$env <- env
  RDAsforAll[[i]]$scores <- scores
}

#save(RDAsforAll, file = "RDAsforAll.RData")
load('RDAsforAll.Rdata')

plotsforAll <- list()
#rdaVecslist <- list()
for(i in systemSeason_list$systemSeason){
 print(i)
  mult <- 1
  mult_vec <- 1
  scores_sites <- data.frame(RDAsforAll[[i]]$scores$sites)
  vecs <- data.frame(RDAsforAll[[i]]$env$vectors$arrows)
  #rdaVecslist[[i]] <- vecs
  middles <- data.frame(RDAsforAll[[i]]$scores$centroids)
  plotSub <- SXS_full %>% #filter out systemSeason of interest just as in rda loop (df_env)
    filter(systemSeason %in% i) %>%
    select(seasonYear) %>%
    cbind(scores_sites) %>%
    group_by(seasonYear) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
  # grps <- data.frame(rownames(middles) %>%
  #   str_remove_all("seasonYear")) %>%
  #   rename(seasonYear = rownames.middles......str_remove_all..seasonYear..)
  
  plotsforAll[[i]] <- ggplot(plotSub) +
    geom_vline(xintercept = 0,
               colour     = "grey70",
               size       = 1) +
    geom_hline(yintercept = 0,
               colour     = "grey70",
               size       = 1) +
    scale_x_continuous(limits       = symmetric_range,
                       breaks       = c(-6,-1,0,1,6),
                       minor_breaks = NULL) +
    scale_y_continuous(limits       = symmetric_range,
                       breaks       = c(-6,-1,0,1,6),
                       minor_breaks = NULL) +
    labs(title = i, 
         x     = paste("CA1 (", round(100*RDAsforAll[[i]]$rda$CCA$eig[1]/RDAsforAll[[i]]$rda$tot.chi,2), "%)", sep = ""),
         y     = paste("CA1 (", round(100*RDAsforAll[[i]]$rda$CCA$eig[2]/RDAsforAll[[i]]$rda$tot.chi,2), "%)", sep = ""))+
    scale_color_manual(values = c("blue","darkgreen"),
                       guide  = guide_none()) +
    # geom_point(data = scores_sites,
    #            aes(x = mult*dbRDA1,
    #                y = mult*dbRDA2,
    #            ),
    #            size   = 1,
    #            stroke = 0.1,
    #            pch    = 21,
    #            colour = "black") +
    geom_point(aes(x = mult*dbRDA1,
                   y = mult*dbRDA2,
                   fill = as.numeric(as.character(seasonYear))),
               size   = 4,
               stroke = 0.1,
               pch    = 21,
               colour = "black") +
    scale_fill_gradient(
      low = "red",
      high = "blue",
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill"
    ) +
    theme_bw() +
    theme(legend.text       = element_text(size=rel(0.8)),
          legend.position   = c(0.055,0.89),
          legend.background = element_blank(),
          legend.key        = element_blank(),
          panel.grid        = element_blank()) +
    geom_segment(data = vecs,
                 aes(x     = 0, 
                     xend  = mult_vec*dbRDA1, 
                     y     = 0, 
                     yend  = mult_vec*dbRDA2, 
                     #color = "blue"
                       ),
                 arrow = arrow(length = unit(.1,"inches")),
                 size  = 1,
                 alpha = 0.6) +
    geom_label_repel(data = vecs,
                     aes(x      = mult_vec*dbRDA1,
                         y      = mult_vec*dbRDA2,
                         label  = rownames(vecs),
                         #colour = "blue",
                         #size   = "blue"
                           ),
                     force         = 0.1,
                     alpha         = 0.8,
                     fontface      = "bold",
                     label.padding = 0.25,
                     box.padding   = 0.1,
                     fill          ="white",
                     label.r       = 0.25,
                     label.size    = 0.25,
                     segment.alpha =0,
                     nudge_x       = ifelse(vecs$dbRDA1>0,0.1,-0.1),
                     nudge_y       = ifelse(vecs$dbRDA2>0,0.2,-0.2)) +
    scale_size_manual(values = c(2.5,3.5),
                      guide  = guide_none()) +
    guides(fill = guide_legend(override.aes = list(size = 2.5),
                               keyheight    = 0.15,
                               keywidth     = 0.2,
                               title = "Year"))
  
}

fullPlot <- wrap_plots(plotsforAll,
                      ncol = 4)

plot(fullPlot)

#ggsave("./Outputs/RDAs/RDAsforAll.tiff", fullPlot, width = 18, height = 9, dpi = 600)

######### BRUTE FORCE PLOTTING ########
# plots$a <- plots$AP_winter +
#   geom_label_repel(data = rdaVecslist$AP_winter,
#                    aes(x      = mult*dbRDA1,
#                        y      = mult*dbRDA2,
#                        label  = rownames(rdaVecslist$AP_winter),
#                        #colour = "blue",
#                        #size   = "blue"
#                    ),
#                    force         = 0.1,
#                    alpha         = 0.8,
#                    fontface      = "bold",
#                    label.padding = 0.25,
#                    box.padding   = 0.1,
#                    fill          ="white",
#                    label.r       = 0.25,
#                    label.size    = 0.25,
#                    segment.alpha =0,
#                    nudge_x       = ifelse(rdaVecslist$AP_winter$dbRDA1>0,0.1,-0.1),
#                    nudge_y       = ifelse(rdaVecslist$AP_winter$dbRDA2>0,0.2,-0.2))
# 
# plots$b <- plots$AP_summer +
#   geom_label_repel(data = rdaVecslist$AP_summer,
#                    aes(x      = mult*dbRDA1,
#                        y      = mult*dbRDA2,
#                        label  = rownames(rdaVecslist$AP_summer),
#                        #colour = "blue",
#                        #size   = "blue"
#                    ),
#                    force         = 0.1,
#                    alpha         = 0.8,
#                    fontface      = "bold",
#                    label.padding = 0.25,
#                    box.padding   = 0.1,
#                    fill          ="white",
#                    label.r       = 0.25,
#                    label.size    = 0.25,
#                    segment.alpha =0,
#                    nudge_x       = ifelse(rdaVecslist$AP_summer$dbRDA1>0,0.1,-0.1),
#                    nudge_y       = ifelse(rdaVecslist$AP_summer$dbRDA2>0,0.2,-0.2))
# 
# plots$c <- plots$CK_winter +
#   geom_label_repel(data = rdaVecslist$CK_winter,
#                    aes(x      = mult*dbRDA1,
#                        y      = mult*dbRDA2,
#                        label  = rownames(rdaVecslist$CK_winter),
#                        #colour = "blue",
#                        #size   = "blue"
#                    ),
#                    force         = 0.1,
#                    alpha         = 0.8,
#                    fontface      = "bold",
#                    label.padding = 0.25,
#                    box.padding   = 0.1,
#                    fill          ="white",
#                    label.r       = 0.25,
#                    label.size    = 0.25,
#                    segment.alpha =0,
#                    nudge_x       = ifelse(rdaVecslist$CK_winter$dbRDA1>0,0.1,-0.1),
#                    nudge_y       = ifelse(rdaVecslist$CK_winter$dbRDA2>0,0.2,-0.2))
# 
# plots$d <- plots$CK_summer +
#   geom_label_repel(data = rdaVecslist$CK_summer,
#                    aes(x      = mult*dbRDA1,
#                        y      = mult*dbRDA2,
#                        label  = rownames(rdaVecslist$CK_summer),
#                        #colour = "blue",
#                        #size   = "blue"
#                    ),
#                    force         = 0.1,
#                    alpha         = 0.8,
#                    fontface      = "bold",
#                    label.padding = 0.25,
#                    box.padding   = 0.1,
#                    fill          ="white",
#                    label.r       = 0.25,
#                    label.size    = 0.25,
#                    segment.alpha =0,
#                    nudge_x       = ifelse(rdaVecslist$CK_summer$dbRDA1>0,0.1,-0.1),
#                    nudge_y       = ifelse(rdaVecslist$CK_summer$dbRDA2>0,0.2,-0.2))
# 
# plots$e <- plots$TB_winter +
#   geom_label_repel(data = rdaVecslist$TB_winter,
#                    aes(x      = mult*dbRDA1,
#                        y      = mult*dbRDA2,
#                        label  = rownames(rdaVecslist$TB_winter),
#                        #colour = "blue",
#                        #size   = "blue"
#                    ),
#                    force         = 0.1,
#                    alpha         = 0.8,
#                    fontface      = "bold",
#                    label.padding = 0.25,
#                    box.padding   = 0.1,
#                    fill          ="white",
#                    label.r       = 0.25,
#                    label.size    = 0.25,
#                    segment.alpha =0,
#                    nudge_x       = ifelse(rdaVecslist$TB_winter$dbRDA1>0,0.1,-0.1),
#                    nudge_y       = ifelse(rdaVecslist$TB_winter$dbRDA2>0,0.2,-0.2))
# 
# plots$f <- plots$TB_summer +
#   geom_label_repel(data = rdaVecslist$TB_summer,
#                    aes(x      = mult*dbRDA1,
#                        y      = mult*dbRDA2,
#                        label  = rownames(rdaVecslist$TB_summer),
#                        #colour = "blue",
#                        #size   = "blue"
#                    ),
#                    force         = 0.1,
#                    alpha         = 0.8,
#                    fontface      = "bold",
#                    label.padding = 0.25,
#                    box.padding   = 0.1,
#                    fill          ="white",
#                    label.r       = 0.25,
#                    label.size    = 0.25,
#                    segment.alpha =0,
#                    nudge_x       = ifelse(rdaVecslist$TB_summer$dbRDA1>0,0.1,-0.1),
#                    nudge_y       = ifelse(rdaVecslist$TB_summer$dbRDA2>0,0.2,-0.2))
# 
# plots$g <- plots$CH_winter +
#   geom_label_repel(data = rdaVecslist$CH_winter,
#                    aes(x      = mult*dbRDA1,
#                        y      = mult*dbRDA2,
#                        label  = rownames(rdaVecslist$CH_winter),
#                        #colour = "blue",
#                        #size   = "blue"
#                    ),
#                    force         = 0.1,
#                    alpha         = 0.8,
#                    fontface      = "bold",
#                    label.padding = 0.25,
#                    box.padding   = 0.1,
#                    fill          ="white",
#                    label.r       = 0.25,
#                    label.size    = 0.25,
#                    segment.alpha =0,
#                    nudge_x       = ifelse(rdaVecslist$CH_winter$dbRDA1>0,0.1,-0.1),
#                    nudge_y       = ifelse(rdaVecslist$CH_winter$dbRDA2>0,0.2,-0.2))
# 
# plots$h <- plots$CH_summer +
#   geom_label_repel(data = rdaVecslist$CH_summer,
#                    aes(x      = mult*dbRDA1,
#                        y      = mult*dbRDA2,
#                        label  = rownames(rdaVecslist$CH_summer),
#                        #colour = "blue",
#                        #size   = "blue"
#                    ),
#                    force         = 0.1,
#                    alpha         = 0.8,
#                    fontface      = "bold",
#                    label.padding = 0.25,
#                    box.padding   = 0.1,
#                    fill          ="white",
#                    label.r       = 0.25,
#                    label.size    = 0.25,
#                    segment.alpha =0,
#                    nudge_x       = ifelse(rdaVecslist$CH_summer$dbRDA1>0,0.1,-0.1),
#                    nudge_y       = ifelse(rdaVecslist$CH_summer$dbRDA2>0,0.2,-0.2))
# 
# 
# 
# testing <- wrap_plots(plots$a,
#                       plots$c,
#                       plots$e,
#                       plots$g,
#                       plots$b,
#                       plots$d,
#                       plots$f,
#                       plots$h,
#                       ncol = 4)
# 
# plot(testing)