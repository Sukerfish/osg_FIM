# CAP analyses for each system/season combination

#### data input ######
library(tidyverse)
library(vegan)
library(BiodiversityR)
#library(MASS)
library(colortools)
library(ggrepel)
library(RColorBrewer)
library(egg)
library(patchwork)

#source https://github.com/jonpeake/forage-fish-public-code/blob/main/ForageFigures_CODE.R
cap_pctvars <- function(cap_out){
  # A function that takes the output of the CAPdiscrim function and 
  # produces the percent of among-group variability explained by each
  # canonical axis.
  #
  # Inputs: 
  # cap_out = 'CAPDiscrim' function list output
  # 
  # Output: 
  # pct_var = A numerical array with dimensions 1 x p, where p is the
  #           number of canonical axes, corresponding to the percent
  #           of among-group variability explained by each axis.
  #         
  
  # Extract eigenvalues from manova sublist
  eig     <- cap_out[["manova"]][["Eigenvalues"]]
  
  # Convert eigenvalues to analog of canonical variability
  vars    <- eig*cap_out$tot*(cap_out$varm/100)/(eig+1)
  
  # Calculate percent of variability explained by each axis
  pct_var <- t(100*vars/sum(vars))
  
  # Trim rows not corresponding to canonical axes
  pct_var <- pct_var[1:ncol(cap_out$x),]
  return(pct_var)
}

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

CAPsforAll <- list()
for(i in systemSeason_list$systemSeason){
  print(i) #watch progress through list
  
  df <- SXS_full %>% #filter out systemSeason of interest
    filter(systemSeason %in% i)
  
  df_spe <- df %>% #pull out taxa only
    subset(select = -c(systemSeason, seasonYear, Reference, systemZone, BottomVegCover, Temperature)) %>%
    select(which(!colSums(., na.rm=TRUE) %in% 0))
  
  df_env <- data.frame(df %>% #pull out environmental variables
    subset(select = c(systemSeason, seasonYear, BottomVegCover, Temperature)) %>%
    mutate(contYear = as.numeric(as.character(seasonYear))))
  
  bf <- (vegdist(df_spe))^0.5 #Bray-Curtis w/ sqrt to reduce negative eigenvalues
  
  CAP <- CAPdiscrim(bf ~ seasonYear,
                    data = df_env,
                    mmax = 200,
                    #add = "lingoes",
                    #parallel = 6,
                    #method="bray", 
                    #permutations = 999
                    )
  
  CAPsforAll[[i]] <- CAP
  CAPsforAll[[i]] <- add.spec.scores(CAPsforAll[[i]], df_spe)

  #cbPalette1    <- brewer.pal(8, "Paired")[c(2, 3, 6, 7)]
  #cbPalette2    <- brewer.pal(4, "Dark2")
}

#save(CAPsforAll, file = "CAPsforAll.RData")
load('CAPSforAll.Rdata')

plots <- list()
sppvecs_list <- list()
for(i in systemSeason_list$systemSeason){
  print(i) #watch progress through list
  var_CA <- round(cap_pctvars(CAPsforAll[[i]]), 2)
  #var_CA <- round(100 * CAPsforAll[[i]][["lda.other"]][["svd"]]^2/sum(CAPsforAll[[i]][["lda.other"]][["svd"]]^2), 2)
  site_scores <- data.frame(CAPsforAll[[i]]$x)
  spp_vecs <- data.frame(CAPsforAll[[i]]$cproj)
  sppvecs_list[[i]] <- spp_vecs
  mult <- 6
  plots[[i]] <- ggplot(site_scores) + 
    geom_vline(xintercept = 0,
               colour     = "grey70",
               size       = .25) +
    geom_hline(yintercept = 0,
               colour     = "grey70",
               size       = .25) +
    scale_x_continuous(limits       = symmetric_range((1+mult)*spp_vecs$LD1),
                       breaks       = c(-6,-3,0,3,6),
                       minor_breaks = NULL) +
    scale_y_continuous(limits       = symmetric_range,
                       breaks       = c(-6,-3,0,3,6),
                       minor_breaks = NULL) +
    geom_point(aes(x    = LD1, 
                   y    = LD2, 
                   fill = NULL),
               size   = 1, 
               stroke = 0.1,
               pch    = 21, 
               colour = "black") +
    labs(title = i,
         x     = paste("CA1 (",var_CA[1],"%)",sep=""),
         y     = paste("CA2 (",var_CA[2],"%)",sep=""),
         fill  = NULL) +
    # scale_fill_manual(values = cbPalette1,
    #                   labels = c(unique(as.character(df_env$seasonYear)))) +
    theme_bw() +
    theme(legend.text       = element_text(size=rel(0.8)),
          legend.position   = c(0.055,0.89),
          legend.background = element_blank(),
          legend.key        = element_blank(),
          panel.grid        = element_blank()) +
    geom_segment(data = spp_vecs,
                 aes(x    = 0,
                     xend = mult*LD1,
                     y    = 0,
                     yend = mult*LD2),
                 arrow = arrow(length = unit(0.1,"inches")),
                 color = "grey1",
                 size  = .75,
                 alpha = 0.6)+
    # geom_label_repel(data = spp_vecs,
    #                  aes(x     = mult*LD1 + 0.2*cos(atan(LD2/LD1))*sign(LD1),
    #                      y     = mult*LD2 + 0.2*sin(atan(LD2/LD1))*sign(LD2),
    #                      label = rownames(spp_vecs)),
    #                  segment.alpha = 0,
    #                  size          = 3.25,
    #                  color         = "blue",
    #                  fontface      = "bold",
    #                  fill          = "white",
    #                  alpha         = 0.7,
    #                  box.padding   = .25,
    #                  lineheight    = 0.4,
    #                  label.size    = 0.25,
    #                  nudge_x       = ifelse(spp_vecs$LD1>0,0.05,-0.05),
    #                  nudge_y       = ifelse(spp_vecs$LD2>-0.02,0.05,-0.05),
    #                  force         = 0.3) +
    guides(fill = guide_legend(override.aes = list(size = 2.5),
                               keyheight    = 0.15,
                               keywidth     = 0.2))
  
  #ggsave(paste("./Outputs/CAPs/CAP_", i, ".tiff", sep = ""), CAPsforAll[[i]]$plots, width = 8, height = 8, dpi = 400)
}


######### BRUTE FORCE PLOTTING ########
plots$a <- plots$AP_winter +
  geom_label_repel(data = sppvecs_list$AP_winter,
                   aes(x     = mult*LD1 + 0.2*cos(atan(LD2/LD1))*sign(LD1),
                       y     = mult*LD2 + 0.2*sin(atan(LD2/LD1))*sign(LD2),
                       label = rownames(sppvecs_list$AP_winter)),
                   segment.alpha = 0,
                   size          = 3.25,
                   color         = "blue",
                   fontface      = "bold",
                   fill          = "white",
                   alpha         = 0.7,
                   box.padding   = .25,
                   lineheight    = 0.4,
                   label.size    = 0.25,
                   nudge_x       = ifelse(sppvecs_list$AP_winter$LD1>0,0.05,-0.05),
                   nudge_y       = ifelse(sppvecs_list$AP_winter$LD2>-0.02,0.05,-0.05),
                   force         = 0.3)

plots$b <- plots$AP_summer +
  geom_label_repel(data = sppvecs_list$AP_summer,
                   aes(x     = mult*LD1 + 0.2*cos(atan(LD2/LD1))*sign(LD1),
                       y     = mult*LD2 + 0.2*sin(atan(LD2/LD1))*sign(LD2),
                       label = rownames(sppvecs_list$AP_summer)),
                   segment.alpha = 0,
                   size          = 3.25,
                   color         = "blue",
                   fontface      = "bold",
                   fill          = "white",
                   alpha         = 0.7,
                   box.padding   = .25,
                   lineheight    = 0.4,
                   label.size    = 0.25,
                   nudge_x       = ifelse(sppvecs_list$AP_summer$LD1>0,0.05,-0.05),
                   nudge_y       = ifelse(sppvecs_list$AP_summer$LD2>-0.02,0.05,-0.05),
                   force         = 0.3)

plots$c <- plots$CK_winter +
  geom_label_repel(data = sppvecs_list$CK_winter,
                   aes(x     = mult*LD1 + 0.2*cos(atan(LD2/LD1))*sign(LD1),
                       y     = mult*LD2 + 0.2*sin(atan(LD2/LD1))*sign(LD2),
                       label = rownames(sppvecs_list$CK_winter)),
                   segment.alpha = 0,
                   size          = 3.25,
                   color         = "blue",
                   fontface      = "bold",
                   fill          = "white",
                   alpha         = 0.7,
                   box.padding   = .25,
                   lineheight    = 0.4,
                   label.size    = 0.25,
                   nudge_x       = ifelse(sppvecs_list$CK_winter$LD1>0,0.05,-0.05),
                   nudge_y       = ifelse(sppvecs_list$CK_winter$LD2>-0.02,0.05,-0.05),
                   force         = 0.3)

plots$d <- plots$CK_summer +
  geom_label_repel(data = sppvecs_list$CK_summer,
                   aes(x     = mult*LD1 + 0.2*cos(atan(LD2/LD1))*sign(LD1),
                       y     = mult*LD2 + 0.2*sin(atan(LD2/LD1))*sign(LD2),
                       label = rownames(sppvecs_list$CK_summer)),
                   segment.alpha = 0,
                   size          = 3.25,
                   color         = "blue",
                   fontface      = "bold",
                   fill          = "white",
                   alpha         = 0.7,
                   box.padding   = .25,
                   lineheight    = 0.4,
                   label.size    = 0.25,
                   nudge_x       = ifelse(sppvecs_list$CK_summer$LD1>0,0.05,-0.05),
                   nudge_y       = ifelse(sppvecs_list$CK_summer$LD2>-0.02,0.05,-0.05),
                   force         = 0.3)

plots$e <- plots$TB_winter +
  geom_label_repel(data = sppvecs_list$TB_winter,
                   aes(x     = mult*LD1 + 0.2*cos(atan(LD2/LD1))*sign(LD1),
                       y     = mult*LD2 + 0.2*sin(atan(LD2/LD1))*sign(LD2),
                       label = rownames(sppvecs_list$TB_winter)),
                   segment.alpha = 0,
                   size          = 3.25,
                   color         = "blue",
                   fontface      = "bold",
                   fill          = "white",
                   alpha         = 0.7,
                   box.padding   = .25,
                   lineheight    = 0.4,
                   label.size    = 0.25,
                   nudge_x       = ifelse(sppvecs_list$TB_winter$LD1>0,0.05,-0.05),
                   nudge_y       = ifelse(sppvecs_list$TB_winter$LD2>-0.02,0.05,-0.05),
                   force         = 0.3)

plots$f <- plots$TB_summer +
  geom_label_repel(data = sppvecs_list$TB_summer,
                   aes(x     = mult*LD1 + 0.2*cos(atan(LD2/LD1))*sign(LD1),
                       y     = mult*LD2 + 0.2*sin(atan(LD2/LD1))*sign(LD2),
                       label = rownames(sppvecs_list$TB_summer)),
                   segment.alpha = 0,
                   size          = 3.25,
                   color         = "blue",
                   fontface      = "bold",
                   fill          = "white",
                   alpha         = 0.7,
                   box.padding   = .25,
                   lineheight    = 0.4,
                   label.size    = 0.25,
                   nudge_x       = ifelse(sppvecs_list$TB_summer$LD1>0,0.05,-0.05),
                   nudge_y       = ifelse(sppvecs_list$TB_summer$LD2>-0.02,0.05,-0.05),
                   force         = 0.3)

plots$g <- plots$CH_winter +
  geom_label_repel(data = sppvecs_list$CH_winter,
                   aes(x     = mult*LD1 + 0.2*cos(atan(LD2/LD1))*sign(LD1),
                       y     = mult*LD2 + 0.2*sin(atan(LD2/LD1))*sign(LD2),
                       label = rownames(sppvecs_list$CH_winter)),
                   segment.alpha = 0,
                   size          = 3.25,
                   color         = "blue",
                   fontface      = "bold",
                   fill          = "white",
                   alpha         = 0.7,
                   box.padding   = .25,
                   lineheight    = 0.4,
                   label.size    = 0.25,
                   nudge_x       = ifelse(sppvecs_list$CH_winter$LD1>0,0.05,-0.05),
                   nudge_y       = ifelse(sppvecs_list$CH_winter$LD2>-0.02,0.05,-0.05),
                   force         = 0.3)

plots$h <- plots$CH_summer +
  geom_label_repel(data = sppvecs_list$CH_summer,
                   aes(x     = mult*LD1 + 0.2*cos(atan(LD2/LD1))*sign(LD1),
                       y     = mult*LD2 + 0.2*sin(atan(LD2/LD1))*sign(LD2),
                       label = rownames(sppvecs_list$CH_summer)),
                   segment.alpha = 0,
                   size          = 3.25,
                   color         = "blue",
                   fontface      = "bold",
                   fill          = "white",
                   alpha         = 0.7,
                   box.padding   = .25,
                   lineheight    = 0.4,
                   label.size    = 0.25,
                   nudge_x       = ifelse(sppvecs_list$CH_summer$LD1>0,0.05,-0.05),
                   nudge_y       = ifelse(sppvecs_list$CH_summer$LD2>-0.02,0.05,-0.05),
                   force         = 0.3)



testing <- wrap_plots(plots$a,
                      plots$c,
                      plots$e,
                      plots$g,
                      plots$b,
                      plots$d,
                      plots$f,
                      plots$h,
                      ncol = 4)

plot(testing)

### centroids ####
for(i in systemSeason_list$systemSeason){
df <- SXS_full %>% #filter out systemSeason of interest
  filter(systemSeason %in% i)

df_spe <- df %>% #pull out taxa only
  subset(select = -c(systemSeason, seasonYear, Reference, systemZone, BottomVegCover, Temperature)) %>%
  select(which(!colSums(., na.rm=TRUE) %in% 0))

df_env <- data.frame(df %>% #pull out environmental variables
                       subset(select = c(systemSeason, seasonYear, BottomVegCover, Temperature)) %>%
                       mutate(contYear = as.numeric(as.character(seasonYear))))

CAPsforAll[[i]]$centroids <- centroids.long(sites.long(ordiplot(CAPsforAll[[i]])),
                                            df_env$seasonYear,
                                            centroids.only = TRUE)
}


##### centroids with sppvec plots #######
for(i in systemSeason_list$systemSeason){
  print(i) #watch progress through list
  var_CA <- round(cap_pctvars(CAPsforAll[[i]]), 2)
  #var_CA <- round(100 * CAPsforAll[[i]][["lda.other"]][["svd"]]^2/sum(CAPsforAll[[i]][["lda.other"]][["svd"]]^2), 2)
  site_scores <- data.frame(CAPsforAll[[i]]$x)
  middles <- data.frame(CAPsforAll[[i]]$centroids)
  spp_vecs <- data.frame(CAPsforAll[[i]]$cproj) %>%
                           arrange(desc(abs(LD1))) %>%
                           slice_head(n = 10)
  sppvecs_list[[i]] <- spp_vecs
  mult <- 6
  plots[[i]] <- ggplot(middles) + 
    geom_vline(xintercept = 0,
               colour     = "grey70",
               size       = .25) +
    geom_hline(yintercept = 0,
               colour     = "grey70",
               size       = .25) +
    scale_x_continuous(limits       = symmetric_range((1+mult)*spp_vecs$LD1),
                       breaks       = c(-6,-3,0,3,6),
                       minor_breaks = NULL) +
    scale_y_continuous(limits       = symmetric_range,
                       breaks       = c(-6,-3,0,3,6),
                       minor_breaks = NULL) +
    geom_point(aes(x    = axis1c, 
                   y    = axis2c, 
                   fill = as.numeric(as.character(df_env.seasonYear))),
               size   = 4, 
               stroke = 0.1,
               pch    = 21, 
               colour = "black") +
    labs(title = i,
         x     = paste("CA1 (",var_CA[1],"%)",sep=""),
         y     = paste("CA2 (",var_CA[2],"%)",sep=""),
         fill  = NULL) +
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
    geom_segment(data = spp_vecs,
                 aes(x    = 0,
                     xend = mult*LD1,
                     y    = 0,
                     yend = mult*LD2),
                 arrow = arrow(length = unit(0.1,"inches")),
                 color = "grey1",
                 size  = .75,
                 alpha = 0.6)+
  # geom_label_repel(data = spp_vecs,
  #                  aes(x     = mult*LD1 + 0.2*cos(atan(LD2/LD1))*sign(LD1),
  #                      y     = mult*LD2 + 0.2*sin(atan(LD2/LD1))*sign(LD2),
  #                      label = rownames(spp_vecs)),
  #                  segment.alpha = 0,
  #                  size          = 3.25,
  #                  color         = "blue",
  #                  fontface      = "bold",
  #                  fill          = "white",
  #                  alpha         = 0.7,
  #                  box.padding   = .25,
  #                  lineheight    = 0.4,
  #                  label.size    = 0.25,
  #                  nudge_x       = ifelse(spp_vecs$LD1>0,0.05,-0.05),
  #                  nudge_y       = ifelse(spp_vecs$LD2>-0.02,0.05,-0.05),
  #                  force         = 0.3) +
  guides(fill = guide_legend(override.aes = list(size = 2.5),
                             keyheight    = 0.15,
                             keywidth     = 0.2))
  
  #ggsave(paste("./Outputs/CAPs/CAP_", i, ".tiff", sep = ""), CAPsforAll[[i]]$plots, width = 8, height = 8, dpi = 400)
}