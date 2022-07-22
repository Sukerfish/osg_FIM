# Bray-Curtis similarity pairwise among all years with respect to first year of data
# Using average annual density of each taxa within each estuary and season
# 
# Model change in abundance through time for each taxa
# 
# Calculate coefficient of the slope of each
# 
# Potential coefficient distributions
#### data input ######

library(tidyverse)

load('TidyGearCode20.Rdata')

CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, systemZone, Scientificname, N2)) %>%
  #filter(season == "summer" | season == "winter")
  filter(season == "winter")
  #filter(season == "summer")

###### main ######
SiteXSpeciesFull <- CleanHauls %>%
  mutate(N2 = N2^.25) %>% #fourth-root transform
  group_by(Reference) %>%
  spread(Scientificname,N2) %>%
  ungroup() %>%
  subset(select = -c(Reference, season, systemZone)) %>%
  replace(is.na(.), 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  group_by(system, seasonYear) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

# TBFull <- SiteXSpeciesFull %>%
#   filter(system == "TB") %>%
#   subset(select = -c(system)) %>%
#   column_to_rownames(var = "seasonYear")

systemLogic <- unique(SiteXSpeciesFull$system)

out <- list()

library(vegan)

for (i in 1:length(systemLogic)){
  tempDF <- SiteXSpeciesFull %>%
  filter(system == as.character(systemLogic[i])) %>%
  subset(select = -c(system)) %>%
  column_to_rownames(var = "seasonYear")
  
  distmat <- as.matrix(vegdist(tempDF, method ="bray"))
  
  name <- as.character(systemLogic[i])
  
  out[[name]] <- distmat[,1]
}

yearLogic <- SiteXSpeciesFull %>%
  subset(select = c(system, seasonYear)) %>%
  unique()

lessgo <- stack(out) %>%
  rename(commsim = values) %>%
  rename(system = ind) %>%
  cbind(seasonYear = yearLogic$seasonYear) %>%
  filter(commsim > 0) %>%
  mutate(season = "winter")


CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, systemZone, Scientificname, N2)) %>%
  #filter(season == "summer" | season == "winter")
  #filter(season == "winter")
  filter(season == "summer")

SiteXSpeciesFull <- CleanHauls %>%
  mutate(N2 = N2^.25) %>% #fourth-root transform
  group_by(Reference) %>%
  spread(Scientificname,N2) %>%
  ungroup() %>%
  subset(select = -c(Reference, season, systemZone)) %>%
  replace(is.na(.), 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  group_by(system, seasonYear) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

# TBFull <- SiteXSpeciesFull %>%
#   filter(system == "TB") %>%
#   subset(select = -c(system)) %>%
#   column_to_rownames(var = "seasonYear")

systemLogic <- unique(SiteXSpeciesFull$system)

out <- list()

for (i in 1:length(systemLogic)){
  tempDF <- SiteXSpeciesFull %>%
    filter(system == as.character(systemLogic[i])) %>%
    subset(select = -c(system)) %>%
    column_to_rownames(var = "seasonYear")
  
  distmat <- as.matrix(vegdist(tempDF, method ="bray"))
  
  name <- as.character(systemLogic[i])
  
  out[[name]] <- distmat[,1]
}

yearLogic <- SiteXSpeciesFull %>%
  subset(select = c(system, seasonYear)) %>%
  unique()

temp <- stack(out) %>%
  rename(commsim = values) %>%
  rename(system = ind) %>%
  cbind(seasonYear = yearLogic$seasonYear) %>%
  filter(commsim > 0) %>%
  mutate(season = "summer")

lessgo <- lessgo %>%
  rbind(temp)

# TBYear <- TBFull %>%
#   subset(select = c(seasonYear))
# 
# TBSpecies <- TBFull %>%
#   subset(select = -c(seasonYear))

# TB_bray = vegdist(TBFull, method='bray')
# distmat <- as.matrix(TB_bray)
# test <- distmat[,1]
# test2 <- test[-1]
# 
# test3 <- TBYear %>%
#   filter(seasonYear != "1998") %>%
#   mutate(commsim = test2)

library(ggplot2)

ggplot(data=lessgo,
       aes(seasonYear, commsim, color = system)) +
  theme_bw() + 
  ggtitle("Bray-Curtis community dissimilarity over time") +
  xlab("Year") +
  scale_x_continuous(breaks = seq(1998, 2020, 2)) +
  ylab("Community dissimilarity value") +
  geom_line() +
  facet_wrap(as.factor(lessgo$season)) +
  theme(strip.text.x = element_text(size = 20)) +
  theme(strip.text.y = element_text(size = 20))

##### adding zone ######


SiteXSpeciesFull <- CleanHauls %>%
  mutate(N2 = N2^.25) %>% #fourth-root transform
  group_by(Reference) %>%
  spread(Scientificname,N2) %>%
  ungroup() %>%
  subset(select = -c(Reference, season)) %>%
  replace(is.na(.), 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  group_by(system, seasonYear, systemZone) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

# TBFull <- SiteXSpeciesFull %>%
#   filter(system == "TB") %>%
#   subset(select = -c(system)) %>%
#   column_to_rownames(var = "seasonYear")

systemLogic <- unique(SiteXSpeciesFull$systemZone)

out <- list()

library(vegan)

for (i in 1:length(systemLogic)){
  tempDF <- SiteXSpeciesFull %>%
    filter(systemZone == as.character(systemLogic[i])) %>%
    subset(select = -c(system, systemZone)) %>%
    column_to_rownames(var = "seasonYear")
  
  distmat <- as.matrix(vegdist(tempDF, method ="bray"))
  
  name <- as.character(systemLogic[i])
  
  out[[name]] <- distmat[,1]
}

yearLogic <- SiteXSpeciesFull %>%
  subset(select = c(system, seasonYear, systemZone)) %>%
  unique() %>%
  arrange(system, systemZone)
  # arrange(factor(system, levels = c("AP", "CK", "CH", "TB")), systemZone)

lessgo <- stack(out) %>%
  rename(commsim = values) %>%
  rename(systemZone = ind) %>%
  cbind(seasonYear = yearLogic$seasonYear) %>%
  cbind(system = yearLogic$system) %>%
  filter(commsim > 0) %>%
  group_by(system, seasonYear) %>%
  summarise(commsim = mean(commsim)) %>%
  mutate(season = "winter")


CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, systemZone, Scientificname, N2)) %>%
  #filter(season == "summer" | season == "winter")
  #filter(season == "winter")
  filter(season == "summer")

SiteXSpeciesFull <- CleanHauls %>%
  mutate(N2 = N2^.25) %>% #fourth-root transform
  group_by(Reference) %>%
  spread(Scientificname,N2) %>%
  ungroup() %>%
  subset(select = -c(Reference, season)) %>%
  replace(is.na(.), 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  group_by(system, seasonYear, systemZone) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

# TBFull <- SiteXSpeciesFull %>%
#   filter(system == "TB") %>%
#   subset(select = -c(system)) %>%
#   column_to_rownames(var = "seasonYear")


systemLogic <- unique(SiteXSpeciesFull$systemZone)

out <- list()

for (i in 1:length(systemLogic)){
  tempDF <- SiteXSpeciesFull %>%
    filter(systemZone == as.character(systemLogic[i])) %>%
    subset(select = -c(system, systemZone)) %>%
    column_to_rownames(var = "seasonYear")
  
  distmat <- as.matrix(vegdist(tempDF, method ="bray"))
  
  name <- as.character(systemLogic[i])
  
  out[[name]] <- distmat[,1]
}

yearLogic <- SiteXSpeciesFull %>%
  subset(select = c(system, seasonYear, systemZone)) %>%
  unique() %>%
  arrange(system, systemZone)
# arrange(factor(system, levels = c("AP", "CK", "CH", "TB")), systemZone)

temp <- stack(out) %>%
  rename(commsim = values) %>%
  rename(systemZone = ind) %>%
  cbind(seasonYear = yearLogic$seasonYear) %>%
  cbind(system = yearLogic$system) %>%
  filter(commsim > 0) %>%
  group_by(system, seasonYear) %>%
  summarise(commsim = mean(commsim)) %>%
  mutate(season = "summer")

fullLessgo <- lessgo %>%
  rbind(temp) %>%
  arrange(factor(system, levels = c("AP", "CK", "CH", "TB")))



ggplot(data=fullLessgo,
       aes(seasonYear, commsim, color = system)) +
  theme_bw() + 
  ggtitle("Bray-Curtis community dissimilarity values over time (averaged by zone)") +
  xlab("Year") +
  scale_x_continuous(breaks = seq(1998, 2020, 2)) +
  ylab("Community dissimilarity value") +
  geom_line() +
  facet_wrap(as.factor(fullLessgo$season)) +
  theme(strip.text.x = element_text(size = 20)) +
  theme(strip.text.y = element_text(size = 20))