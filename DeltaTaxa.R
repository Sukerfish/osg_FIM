# Heatmap of zscored abundance data
#
#
#### data input ######

library(viridis)
library(tidyverse)
library(ggplot2)

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
  filter(LTmean > 0) %>% #remove taxa entirely absent from each system - if LTmean = 0, taxa never observed
  filter(Scientificname != "No fish") %>%
  #filter(Scientificname == "Leiostomus xanthurus") %>%
  filter(system == "AP") %>%
  filter(season == "winter")

#### plot heatmaps ####

#calculate zscore limits
maxZS <- max(abs(YearXSpeciesZ$zscore), na.rm = TRUE)

ggplot(YearXSpeciesZ, aes(seasonYear, Scientificname, fill= zscore)) + 
  scale_fill_gradientn(colours=c("blue","white", "red"), 
                       na.value = "grey98",
                       limits = c(maxZS*-1, maxZS),
                       ) +
  geom_tile()

ggplot(centered, aes(x=seasonYear, 
                     y=zscore)) + 
  geom_line() +
  facet_wrap(vars(system))

#### Dornelas style ####
#
#
# Linear regression to population abundances
#   ignoring when species was absent (pre/post)
# Sqrt first
# Then scale for mean 0 and stdev 1
# Fit OLS regression through trans data and calc slope/stat sig


# binaryAbs <- rawAbs %>%
#   mutate(logic = if_else(N2 > 0, 1, 0)) %>%
#   group_by(system, season, seasonYear, Scientificname) %>%
#   summarise(avg = mean(logic))

# SiteXSpeciesFull <- CleanHauls %>%
#   #mutate(N2 = N2^.5) %>% #square root transform
#   mutate(logic = if_else(N2 > 0, 1, 0)) %>%
#   #group_by(Reference) %>%
# pivot_wider(id_cols = Reference:systemZone,
#             #id_expand = TRUE,
#             names_from = Scientificname,
#             values_from = logic,
#             values_fill = 0) %>% #replace all NA values with 0s, i.e. counting as true zero
#   ungroup() %>%
#   subset(select = -c(Reference, systemZone)) %>%
#   group_by(system, season, seasonYear) %>%
#   summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
#   group_by(system, season, seasonYear) %>%
#   mutate(across(everything(), ~ if_else(.x > 0, 1, 0))) %>%
# filter(system == "TB") %>%
# filter(season == "summer")

library(tseries)

binaryAbs <- YearXSpeciesZ %>%
  mutate(logic = if_else(avg > 0, 1, 0)) %>%
  pivot_wider(id_cols = system:seasonYear,
              names_from = Scientificname,
              values_from = logic,
              values_fill = 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  filter(system == "CK") %>%
  filter(season == "summer")

seasonLogic <- TidyRefsList %>%
  select(system:season) %>%
  distinct() %>%
  unite(seasonLogic, 
        sep = "_")

runs.p <- runs.test(as.factor(binaryAbs$`Hypsoblennius hentz`),
          alternative = "less")

N.turnovers <- function (vec=rbinom(50,1,0.5)) {
  
  z <- diff(vec)   # get difference lag one
  col <- sum(z==1) # count 0 -> 1
  ext <- sum(z==-1)# count 1 -> 0
  
  result <- c(col,ext)
  names(result) <- c("col","ext")
  return(result)
  
}

WLEC<-data.frame(ID=0,Species=0,runs.pvalue=0,col=0,ext=0,slope=0,slope.pvalue=0,Ninit=0,meanN=0)

binaryAbsTaxa <- binaryAbs %>%
  ungroup() %>%
  select(!c(system:seasonYear))

for(i in 1:ncol(binaryAbsTaxa)){
  df <- binaryAbs
  
  z <- binaryAbs[[i+3]]
  try(runs.p <- runs.test(as.factor(z),alternative="less"), silent = TRUE)
  turnover <- N.turnovers(z)
  
  #WLEC[i,1] <- 
  WLEC[i,2] <- colnames(binaryAbs[i+3])
  WLEC[i,3] <- runs.p$p.value
  WLEC[i,4] <- turnover[1]
  WLEC[i,5] <- turnover[2]
  
}

WLEC <- WLEC %>%
  mutate(test = ifelse(runs.pvalue < 0.05, 1, 0)) %>%
  filter(col > 0 & test == 1)





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