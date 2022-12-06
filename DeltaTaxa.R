# Heatmap of zscored abundance data
#
#
#### data input ######

library(viridis)
library(tidyverse)
library(ggplot2)
library(writexl)

load('TidyGearCode20.Rdata')

#### Z scored abundance ####

#grab basic haul data and counts
CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, systemZone, Scientificname, N2))

#full site by species matrix
SiteXSpeciesFull <- CleanHauls %>%
  pivot_wider(id_cols = Reference:systemZone,
              names_from = Scientificname,
              values_from = N2,
              values_fill = 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  ungroup()

#collapse full site x species matrix to long term means
LTxSpeciesMeans <- SiteXSpeciesFull %>%
  subset(select = -c(Reference, systemZone, seasonYear)) %>%
  group_by(system, season) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

#make matching standard deviation matrix to long term means
LTxSpeciesStdev <- SiteXSpeciesFull %>%
  subset(select = -c(Reference, systemZone, seasonYear)) %>%
  group_by(system, season) %>%
  summarise(across(everything(), ~ sd(.x, na.rm = TRUE)))

#Z score conversion process
YearXSpeciesZ <- SiteXSpeciesFull %>%
  subset(select = -c(Reference, systemZone)) %>%
  group_by(system, season, seasonYear) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>% #collapse to annual means 
  ungroup() %>%
  pivot_longer(cols = !c(system:seasonYear),
              names_to = "Scientificname",
              values_to = "avg") %>%
              #expand back out to long form for leftjoins with LT mean and LT stdev
  ungroup() %>%
  left_join(pivot_longer(data = LTxSpeciesMeans, #LT mean column added
                         cols = !c(system:season),
                         names_to = "Scientificname",
                         values_to = "LTmean")) %>%
  left_join(pivot_longer(data = LTxSpeciesStdev, #LT stdev column added
                         cols = !c(system:season),
                         names_to = "Scientificname",
                         values_to = "stdev")) %>%
  group_by(system, season, seasonYear) %>%
  mutate(zscore = ((avg - LTmean)/stdev)) %>% #calculate zscores using annual means
  ungroup() %>%
  filter(LTmean > 0) %>% #remove taxa entirely absent from each system - if LTmean = 0, taxa never observed
  filter(Scientificname != "No fish") 
  # #filter(Scientificname == "Leiostomus xanthurus") %>%
  # filter(system == "AP") %>%
  # filter(season == "winter")

#### plot heatmaps ####

#calculate zscore limits
# maxZS <- max(abs(YearXSpeciesZ$zscore), na.rm = TRUE)
# 
# ggplot(YearXSpeciesZ, aes(seasonYear, Scientificname, fill= zscore)) + 
#   scale_fill_gradientn(colours=c("blue","white", "red"), 
#                        na.value = "grey98",
#                        limits = c(maxZS*-1, maxZS),
#                        ) +
#   geom_tile()

# ggplot(centered, aes(x=seasonYear, 
#                      y=zscore)) + 
#   geom_line() +
#   facet_wrap(vars(system))

#### Dornelas style ####
#
#
# Linear regression to population abundances
#   ignoring when species was absent (pre/post)
# Sqrt first
# Then scale for mean 0 and stdev 1
# Fit OLS regression through trans data and calc slope/stat sig


library(tseries)

N.turnovers <- function (vec=rbinom(50,1,0.5)) {
  
  z <- diff(vec)   # get difference lag one
  col <- sum(z==1) # count 0 -> 1
  ext <- sum(z==-1)# count 1 -> 0
  
  result <- c(col,ext)
  names(result) <- c("col","ext")
  return(result)
  
}

# binaryAbs <- YearXSpeciesZ %>%
#   mutate(logic = if_else(avg > 0, 1, 0)) %>%
#   pivot_wider(id_cols = system:seasonYear,
#               names_from = Scientificname,
#               values_from = logic,
#               values_fill = 0) %>% #replace all NA values with 0s, i.e. counting as true zero
#   filter(system == "CH") %>%
#   filter(season == "winter")

# seasonLogic <- TidyRefsList %>%
#   select(system:season) %>%
#   distinct() %>%
#   unite(seasonLogic,
#         sep = "_")
#
# runs.p <- runs.test(as.factor(binaryAbs$`Hypsoblennius hentz`),
#           alternative = "less")
# 
# binaryAbsTaxa <- binaryAbs %>%
#   ungroup() %>%
#   select(!c(system:seasonYear))

##### FIRST RUN WITH XFORM DATA FOR SIGNIFICANCE TEST #######

SiteXSpeciesAnnual <- YearXSpeciesZ %>%
  mutate(avgXform = avg^.5) %>% #sqrt transform to account for overly abundant taxa
  pivot_wider(id_cols = system:seasonYear,
              names_from = Scientificname,
              values_from = avgXform,
              values_fill = 0) #replace all NA values with 0s, i.e. counting as true zero

SiteXSpeciesList <- SiteXSpeciesAnnual %>%
  unite("systemSeason",
        system:season,
        sep = "_") %>%
  #mutate(systemSeason = factor(systemSeason, levels = unique(systemSeason))) %>%
  split(.$systemSeason)
  
systemSeason <- SiteXSpeciesAnnual %>%
  select(system:season) %>%
  unite("systemSeason",
        system:season,
        sep = "_") %>%
  distinct()

SlopesForAll <- list() #initialize
for(i in 1:nrow(systemSeason)){
  df <- SiteXSpeciesList[[i]] %>%
    select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) #removes taxa entirely absent from each systemSeason
  zf <- data.frame(Scientificname=0,slope=0,stderr=0,p.value=0) #initialize rolling slope output df
  label <- unique(df$systemSeason) #pull out identifier
  
  for(j in 1:(ncol(df)-2)){
    zf[,1]<-colnames(df[j+2]) #grab column name (i.e. Scientific name) and drop in the row
    zf[,2:4]<-coef(summary(lm(df[[j+2]]~df$seasonYear)))[2,c(1,2,4)] #run lm and extract coef, p-value, SE
    
    #zf[,2:4]<-try(coef(summary(glm(df[[j+2]]~df$seasonYear, family="poisson")))[2,c(1,2,4)], silent = TRUE)
    
    ifelse(j == 1, 
           SlopesForAll[[label]] <- zf, #catch initial "incorrect number of subscripts on matrix" issue
           SlopesForAll[[label]][j,1:4] <- zf) #proceed with all others
  }
}

SlopesForAllDF <- bind_rows(SlopesForAll, .id = "systemSeason") %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_") %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))

###### RERUN WITH RAW AB FOR TRUE SLOPE VALUE (in separate DFs) ######

SiteXSpeciesAnnualRaw <- YearXSpeciesZ %>%
  #mutate(avgXform = avg^.5) %>% #sqrt transform to account for overly abundant taxa
  pivot_wider(id_cols = system:seasonYear,
              names_from = Scientificname,
              values_from = avg,
              values_fill = 0) #replace all NA values with 0s, i.e. counting as true zero

SiteXSpeciesListRaw <- SiteXSpeciesAnnualRaw %>%
  unite("systemSeason",
        system:season,
        sep = "_") %>%
  #mutate(systemSeason = factor(systemSeason, levels = unique(systemSeason))) %>%
  split(.$systemSeason)

systemSeasonRaw <- SiteXSpeciesAnnualRaw %>%
  select(system:season) %>%
  unite("systemSeason",
        system:season,
        sep = "_") %>%
  distinct()

SlopesForAllRaw <- list()
for(i in 1:nrow(systemSeasonRaw)){
  df <- SiteXSpeciesListRaw[[i]] %>%
    select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) #removes taxa entirely absent from each systemSeason
  zf <- data.frame(Scientificname=0,raw_slope=0,stderr=0,p.value=0) #initialize rolling slope output df
  label <- unique(df$systemSeason) #pull out identifier
  
  for(j in 1:(ncol(df)-2)){
    zf[,1]<-colnames(df[j+2]) #grab column name (i.e. Scientific name) and drop in the row
    zf[,2:4]<-coef(summary(lm(df[[j+2]]~df$seasonYear)))[2,c(1,2,4)] #run lm and extract coef, p-value, SE
    
    #zf[,2:4]<-try(coef(summary(glm(df[[j+2]]~df$seasonYear, family="poisson")))[2,c(1,2,4)], silent = TRUE)
    
    ifelse(j == 1, 
           SlopesForAllRaw[[label]] <- zf, #catch initial "incorrect number of subscripts on matrix" issue
           SlopesForAllRaw[[label]][j,1:4] <- zf) #proceed with all others
  }
}

SlopesForAllDFRaw <- bind_rows(SlopesForAllRaw, .id = "systemSeason") %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_") %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))

#add raw slope into original DF
SlopesForAllDF$raw_slope <- SlopesForAllDFRaw$raw_slope 
SlopesForAllDF$raw_stderr <- SlopesForAllDFRaw$stderr

###### plot it up ######
ggplot(SlopesForAllDF, 
       aes(raw_slope))+
  geom_histogram(binwidth = .01) +
  geom_density(adjust = 10, fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  #geom_vline(aes(xintercept = mean), test, linetype="dashed", colour = "blue")+
  geom_vline(xintercept = 0, linetype="dashed") +
  #geom_vline(aes(xintercept = mean+siq), test, colour = "red")+
  #geom_vline(aes(xintercept = mean-siq), test, colour = "red")+
  facet_grid(season~system,
             scales = "free_y") +
  coord_cartesian(xlim = c(-0.05, 0.05)) +
  #xlab("Population Change") +
  #ylab("Number of Taxa") +
  theme(axis.text=element_text(size = 12)) +
  theme(axis.title=element_text(size = 16)) +
  theme(strip.text = element_text(size = 16)) +
  #ggtitle("GLM w/ no Xform") +
  theme(title=element_text(size = 20))

##### significant slopes wrangling ######

# SlopeStats <- SlopesForAllDF %>%
#   group_by(system, season) %>% 
#   summarise(stdev = sd(slope),
#             mean = mean(slope))

sigSlopes <- SlopesForAllDF %>%
  group_by(system, season) %>%
  filter(p.value < 0.05) %>%
  arrange(system, season) %>%
  unite("systemSeason",
        system:season,
        sep = "_")

#make several uniquely ordered plots and bind them...
library(patchwork)

#get unique list for the loop
seasysKey <- sigSlopes %>%
 distinct(systemSeason)

plots <- list() #initialize
#use actual values from the key list
for (i in seasysKey$systemSeason){
  #filter by key values
  plotup <- sigSlopes %>%
    filter(systemSeason == i)
  
  #make list of plots from the filtered plotup DF
  plots[[i]] = ggplot(plotup,
         aes(y = reorder(Scientificname, raw_slope), #order by decreasing slope top to bottom
             x = raw_slope)) + 
    geom_point(aes(color = cut(raw_slope, c(-Inf, 0, Inf))), #cut the slope data for +/- colors
                           size = 3) +
    scale_color_manual(name = "slope", #define the +/- colors
                       values = c("(0, Inf]" = "blue",
                                  "(-Inf,0]" = "red"),
                       labels = c("pos", "neg")) +
    geom_errorbarh(aes(xmin=(raw_slope + (-2*stderr)), xmax=(raw_slope + (2*stderr))),
                   color = "darkgray") + #std err bars
    ggtitle(i) + #dynamic title using the key value
    theme(axis.title.y=element_blank()) +
    coord_cartesian(xlim = c(-.2, 0.2)) +
    geom_vline(xintercept = 0, linetype="dashed")
}

# png(file = "~/osg_FIM/Outputs/sigslopesrawSxS.png",
#     width = 1920, height = 1080)

#wrap them from the list
wrap_plots(plots, ncol = 2, guides = "collect")

# dev.off()

library(rfishbase)

#call FishBase
sigEcoFull <- ecosystem(species_list = unique(c(sigSlopes$Scientificname)))

sigEco <- sigEcoFull %>%
  select(Species, Climate) %>%
  group_by(Species) %>%
  distinct() %>% #parse out unique species/climate combos
  mutate(ClimateClean = str_to_lower(str_squish(Climate))) %>% #clean them up
  select(Species, ClimateClean) %>%
  group_by(Species) %>%
  distinct() %>% #parse out unique combos again
  filter(!is.na(ClimateClean)) %>% #remove nonexistent ones for now
  mutate(presence = 1) %>% #prepare for wide pivot
  pivot_wider(names_from = ClimateClean, values_from = presence) %>%
  #bitfield coding setup (tropical is inherently value 1)
  mutate(subtropical = subtropical*2) %>% 
  mutate(temperate = temperate*4) %>%
  mutate(boreal = NA) %>%
  #sum bitfield rowwise to compress to one data field
  rowwise() %>%
  mutate(tempLogic = sum(subtropical, tropical, temperate, boreal, na.rm = TRUE)) %>%
  #coerce as factor and recode using key
  mutate(tempLogic = as.factor(tempLogic)) %>%
  mutate(tempLogic = recode_factor(tempLogic,
                                   '1' = "tropical",
                                   '2' = "subtropical",
                                   '3' = "tropical/subtropical",
                                   '4' = "temperate",
                                   '5' = "tropical/temperate",
                                   '6' = "subtropical/temperate",
                                   '7' = "tropical/subtropical/temperate")) %>%
  select(Species, tempLogic) %>% #select and rename for leftjoin with other data
  rename(Scientificname = Species)

sigSlopesEco <- sigSlopes %>%
  left_join(sigEco) %>%
  filter(!is.na(tempLogic)) 
 
plotsEco <- list() #initialize
#use actual values from the key list
for (i in seasysKey$systemSeason){
  #filter by key values
  plotup <- sigSlopesEco %>%
    filter(systemSeason == i)
  
  #make list of plots from the filtered plotup DF
  plotsEco[[i]] = ggplot(plotup,
                      aes(y = reorder(Scientificname, raw_slope), #order by decreasing slope top to bottom
                          x = raw_slope,
                          shape = tempLogic)) + 
    geom_point(aes(color = cut(raw_slope, c(-Inf, 0, Inf))), #cut the slope data for +/- colors
               size = 3) +
    scale_color_manual(name = "slope", #define the +/- colors
                       values = c("(0, Inf]" = "blue",
                                  "(-Inf,0]" = "red"),
                       labels = c("pos", "neg")) +
    geom_errorbarh(aes(xmin=(raw_slope + (-2*stderr)), xmax=(raw_slope + (2*stderr))),
                   color = "darkgray") + #std err bars
    scale_shape_manual(
      breaks = c(
        "tropical",
        "subtropical",
        "tropical/subtropical",
        "temperate",
        "tropical/temperate",
        "subtropical/temperate",
        "tropical/subtropical/temperate"
      ),
      values = c(15, 16, 17, 18, 0, 1, 2), #use scale shape identities
      drop = FALSE
    ) +
    ggtitle(i) + #dynamic title using the key value
    theme(axis.title.y=element_blank()) +
    coord_cartesian(xlim = c(-.2, 0.2)) +
    geom_vline(xintercept = 0, linetype="dashed")
}

# png(file = "~/osg_FIM/Outputs/sigslopesrawSxSwithEco.png",
#     width = 1920, height = 1080)

wrap_plots(plotsEco, ncol = 2, guides = "collect")

# dev.off()



######## old stuff #########

# SOI = c("summer",
#         "winter")
# 
# YOI = c("AP",
#         "CK",
#         "TB",
#         "CH")
# 
# SlopesForAll<-data.frame(ID=0,Species=0,runs.pvalue=0,col=0,ext=0,slope=0,slope.pvalue=0,Ninit=0,meanN=0)
# SlopesForAll <- list()
# 
# for(i in length(SOI)){
#   df <- SiteXSpeciesAnnual %>%
#     filter(season == SOI[i])
#   
#   for(j in length(YOI)){
#     zf <- df %>%
#       filter(system == YOI[j]) %>%
#       select(all_of(names(.)[1:3]), where(~ is.numeric(.) && 
#                                             sum(., na.rm = TRUE) > 0))
#     SlopesForAll[[j,2]]   <- paste(colnames(zf[j+3]))
#     #SlopesForAll[[j,2]]   <- colnames(zf[j+3])
#     #SlopesForAll[[j,6]] <- coef(summary(lm(zf[[j+3]]~zf$seasonYear)))[2,c(1)]
#   }
# }
#   
#   #z <- binaryAbs[[i+3]]
#   #runs.p <- list()
#   #try(runs.p <- runs.test(as.factor(z),alternative="less"), silent = TRUE)
#   turnover <- N.turnovers(z)
# 
#   #SlopesForAll[i,1] <- 
#   SlopesForAll[i,2] <- colnames(binaryAbs[i+3])
#   #SlopesForAll[i,3] <- runs.p$p.value
#   SlopesForAll[i,4] <- turnover[1]
#   SlopesForAll[i,5] <- turnover[2]
#   #SlopesForAll[i,6:7]<-coef(summary(glm(SiteXSpeciesAnnual[[i+3]]~SiteXSpeciesAnnual$seasonYear, family="poisson")))[2,c(1,4)]
#   SlopesForAll[i,6:7]<-coef(summary(lm(SiteXSpeciesAnnual[[i+3]]~SiteXSpeciesAnnual$seasonYear)))[2,c(1,4)]
#   
# }

# test <- SlopesForAll %>%
#   mutate(test = ifelse(slope.pvalue < 0.05, 1, 0)) %>%
#   mutate(testslope = ifelse(slope < 0, -1, 1)) %>%
#   filter(test == 1) %>%
#   filter(testslope == -1)
  #filter(col > 0 & test == 1)

SigSlopes <- SlopesForAllDF %>%
  group_by(system, season) %>%
  filter(p.value < 0.05)
  # filter(system == "AP") %>%
  # filter(season == "summer")

#write_xlsx(SigSlopes, "~/osg_FIM/Outputs/SigSlopes.xlsx")

plotup <- YearXSpeciesZ %>%
  group_by(system, season) %>%
  subset(Scientificname %in% SigSlopes$Scientificname) %>%
  left_join(SigSlopes) %>%
  na.exclude %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH")))

ggplot(plotup, aes(x=seasonYear, 
                     y=avg)) + 
  geom_line(aes(color=Scientificname)) +
  facet_grid(season~system,
             scales = "free_y") +
  ggtitle("Sig Slope Linear avg abund") +
  #theme(title=element_text(size = 20)) +
  theme(legend.position="none")

test <- SlopesForAllDF %>%
  mutate(test = ifelse(p.value < 0.05, 1, 0)) %>%
  mutate(testslope = ifelse(slope < 0, -1, 1)) %>%
  filter(test == 1) %>%
  filter(testslope == 1)
#filter(col > 0 & test == 1)

plotup <- YearXSpeciesZ %>%
  subset(Scientificname %in% test$Scientificname) %>%
  filter(system == "AP") %>%
  filter(season == "winter")

ggplot(plotup, aes(x=seasonYear, 
                   y=zscore,
                   color=Scientificname,
                   group=Scientificname)) + 
  ggtitle("Apalach Winter +") +
  geom_line()








idplace<-1
for(id in idsconstant){
  # getting data for relevant studyID
  data<-TS[TS$ID==id,]
  data<- data[data$Abundance>0,]
  groups<-data.frame(as.character(data$Species),as.numeric(data$Year))
  data.mat<- tapply(data$Abundance,groups, FUN=sum)
  # formating data into species by time matrix
  data.mat[is.na(data.mat)]<-0
  #removing species that are always absent
  N.species<-dim(data.mat)[1] #getting number of species
  WLEC[idplace:(idplace+N.species-1),]<-cbind(rep(id,N.species),rownames(data.mat),rep(NA,N.species),rep(NA,N.species),rep(NA,N.species),rep(NA,N.species),rep(NA,N.species),data.mat[,1],apply(data.mat,1,mean))
  
  ############### identifying extinctions and colonizations
  Binary.Data <- data.mat #using previously simulated data matrix
  Binary.Data[Binary.Data > 0] <- 1 # convert to binary
  
  # remove rows in which there were no absences
  Binary.Data <- Binary.Data[which(rowSums(Binary.Data)< ncol(Binary.Data)),]
  
  # create data frame to hold results
  if(!is.matrix(Binary.Data)) Binary.Data <- t(Binary.Data)# HS added
  
  if(dim(Binary.Data)[1]>0){
    # loop through the data
    for (i in 1:nrow(Binary.Data)) {
      # Extract data for a species, conduct runs test, save output
      z <- Binary.Data[i,]
      runs.p <- runs.test(as.factor(z),alternative="less")
      # Count number of colonizations and extinctions, save output
      turnover <- N.turnovers(z)
      WLEC[WLEC$ID==id & WLEC$Species==rownames(Binary.Data)[i],3:5]<-c(runs.p$p.value, turnover[1], turnover[2])
    }
  }
  ##### Winers and Losers
  data.mat<-t(apply(data.mat,1,scale))# scaling by subtracting mean and dividing by sd
  Time<- unique(data$Year) #getting time vector  
  ######### identifying species with positive (Winers) and negative (Losers) trends
  WLEC[idplace:(idplace+N.species-1),6:7]<-t(apply(data.mat,1,get.coeff,x=Time)) # apply to each row of data
  
  idplace<-idplace+N.species
}