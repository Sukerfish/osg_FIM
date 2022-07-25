# diagnostic histograms for data exploration
#
#
#### data input ######

library(tidyverse)
library(ggplot2)

load('TidyGearCode20.Rdata')

#### richness ####
OnlyFish <- TidyBio %>%
  filter(Scientificname != "No fish")

NoFish <- TidyBio %>%
  filter(Scientificname == "No fish") %>%
  mutate(n = 0) %>%
  subset(select = c(Reference, n, system, season, seasonYear, systemZone))

FullRichness <- OnlyFish %>%
  group_by(Reference) %>%
  count() %>%
  inner_join(TidyRefsList) %>%
  ungroup() %>%
  bind_rows(NoFish)

TotalAbundance <- OnlyFish %>%
  group_by(Reference) %>%
  summarise(n = sum(N2)) %>%
  inner_join(TidyRefsList) %>%
  ungroup() %>%
  bind_rows(NoFish) %>%
  mutate(n = n^.25)

# class(FullRichness$system)
# unique(FullRichness$system)
# as.factor(FullRichness$system)
# 
# , levels = c("AP", "CK", "TB", "CH")
# as.ordered), levels = c("AP", "CK", "TB", "CH")

ggplot(data=FullRichness)+
  geom_histogram(binwidth = 1,
                 aes(x=n)) +
  theme_bw() + 
  ggtitle("Richness per haul by estuary") +
  xlab("Number of taxa in each haul") +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  ylab("Frequency") +
  facet_grid(rows = vars(FullRichness$season), 
             cols = vars(FullRichness$system)) +
  theme(
    text=element_text(size=20))
        
#total abundance
ggplot(data=TotalAbundance)+
  geom_histogram(binwidth = 1,
                 aes(x=n)) +
  theme_bw() + 
  ggtitle("Total abundance per haul by estuary") +
  xlab("Number of individuals in each haul (fourth-root transformed)") +
  scale_x_continuous(breaks = seq(0, 100, 3)) +
  ylab("Frequency") +
  facet_grid(rows = vars(TotalAbundance$season), 
             cols = vars(TotalAbundance$system)) +
  theme(
    text=element_text(size=20))

# p + theme(text=element_text(size=20), #change font size of all text
#           axis.text=element_text(size=20), #change font size of axis text
#           axis.title=element_text(size=20), #change font size of axis titles
#           plot.title=element_text(size=20), #change font size of plot title
#           legend.text=element_text(size=20), #change font size of legend text
#           legend.title=element_text(size=20)) #change font size of legend title 
# 
# TBFull <- SiteXSpeciesFull %>%
#   filter(system == "TB") %>%
#   subset(select = -c(system)) %>%
#   column_to_rownames(var = "seasonYear")

CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, systemZone, Scientificname, N2)) %>%
  #filter(season == "summer" | season == "winter")
  filter(season == "winter")
  #filter(season == "summer")

SiteXSpeciesFull <- CleanHauls %>%
  mutate(N2 = N2^.25) %>% #fourth-root transform
  group_by(Reference) %>%
  spread(Scientificname,N2) %>%
  ungroup() %>%
  subset(select = -c(Reference, season)) %>%
  replace(is.na(.), 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  group_by(system, seasonYear, systemZone) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

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