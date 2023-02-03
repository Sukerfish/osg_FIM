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
  subset(select = c(Reference, system, season, seasonYear, BottomVegCover, systemZone, Scientificname, N2))
  #filter(season == "summer" | season == "winter")
  #filter(season == "winter")
  #filter(season == "summer")

##### NMDS/annual PERMANOVA #######
library(vegan)

SXS_annual <- CleanHauls %>%
  mutate(N2 = N2^.25) %>% #fourth-root transform
  group_by(Reference) %>%
  filter(sum(N2)>0) %>% #remove all References with 0 taxa found
  spread(Scientificname,N2) %>%
  ungroup() %>%
  #subset(select = -c(Reference, season, systemZone)) %>%
  replace(is.na(.), 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  group_by(season, system, seasonYear) %>%
  mutate(seasonYear = as.factor(seasonYear)) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

SXS_winter_spe <- SXS_annual %>%
  filter(season == "winter") %>%
  subset(select = -c(season, system, seasonYear, Reference,systemZone, BottomVegCover,season)) %>%
  select(which(!colSums(., na.rm=TRUE) %in% 0))

SXS_winter_env <- SXS_annual %>%
  filter(season == "winter") %>%
  subset(select = c(system,seasonYear,BottomVegCover))

winterAbundNMDS = metaMDS(SXS_winter_spe,
                          distance = "bray",
                          k=2)
SXS_winter_sppfit <- envfit(winterAbundNMDS, SXS_winter_spe, permutations = 999)

plot(winterAbundNMDS, type = "n", las = 1, main = "NMDS winter")
points(winterAbundNMDS, display = "sites")
points(winterAbundNMDS, display = "species", col = "red", pch = 3)
ordiellipse(winterAbundNMDS,
            groups = SXS_winter_env$system,
            kind = "se",
            conf = 0.95,
            display = "sites",
            label=T)

plot(winterAbundNMDS, main = "NMDS winter")
ordihull(winterAbundNMDS,groups=SXS_winter_env$system,draw="polygon",col="grey90",label=T)
plot(SXS_winter_sppfit, p.max = 0.001, col = "black", cex = 0.7)


winterAbundPERM = adonis2(SXS_winter_spe ~ system+seasonYear+BottomVegCover,
       data = SXS_winter_env,
       method="bray",
       permutations=999)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = SXS_winter_spe ~ system + seasonYear + BottomVegCover, data = SXS_winter_env, permutations = 999, method = "bray")
# Df SumOfSqs      R2       F Pr(>F)    
# system          3   6.4588 0.66082 71.0720  0.001 ***
# seasonYear     22   1.4815 0.15157  2.2230  0.001 ***
# BottomVegCover  1   0.0465 0.00475  1.5338  0.158    
# Residual       59   1.7872 0.18286                   
# Total          85   9.7739 1.00000                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



SXS_summer_spe <- SXS_annual %>%
  filter(season == "summer") %>%
  subset(select = -c(season, system, seasonYear, Reference,systemZone, BottomVegCover,season)) %>%
  select(which(!colSums(., na.rm=TRUE) %in% 0))

SXS_summer_env <- SXS_annual %>%
  filter(season == "winter") %>%
  subset(select = c(system,seasonYear,BottomVegCover))

summerAbundNMDS = metaMDS(SXS_summer_spe,
                          distance = "bray",
                          k=2)
SXS_summer_sppfit <- envfit(summerAbundNMDS, SXS_summer_spe, permutations = 999)

plot(summerAbundNMDS, type = "n", las = 1, main = "NMDS summer")
points(summerAbundNMDS, display = "sites")
points(summerAbundNMDS, display = "species", col = "red", pch = 3)
ordiellipse(summerAbundNMDS,
            groups = SXS_summer_env$system,
            kind = "se",
            conf = 0.95,
            display = "sites",
            label=T)

plot(summerAbundNMDS, main = "NMDS summer")
ordihull(summerAbundNMDS,groups=SXS_summer_env$system,draw="polygon",col="grey90",label=T)
#plot(SXS_summer_sppfit, p.max = 0.001, col = "black", cex = 0.7)


summerAbundPERM = adonis2(SXS_summer_spe ~ system+seasonYear+BottomVegCover,
                          data = SXS_summer_env,
                          method="bray",
                          permutations=999)

plot(summerAbundNMDS, type = "n", las = 1, main = "NMDS summer")
points(summerAbundNMDS, display = "sites")
points(summerAbundNMDS, display = "species", col = "red", pch = 3)
ordihull(summerAbundNMDS,groups = SXS_summer_env$system,display = "sites",label=T)



###### PERMANOVAs pairwise ######
SXS_full <- CleanHauls %>%
  mutate(N2 = N2^.25) %>% #fourth-root transform
  group_by(Reference) %>%
  filter(sum(N2)>0) %>% #remove all References with 0 taxa found
  spread(Scientificname,N2) %>%
  ungroup() %>%
  #subset(select = -c(Reference, season, systemZone)) %>%
  replace(is.na(.), 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  group_by(season, system, seasonYear) %>%
  mutate(seasonYear = as.factor(seasonYear))
  #summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

SXSf_winter_spe <- SXS_full %>%
  filter(season == "winter") %>%
  subset(select = -c(season, system, seasonYear, Reference,systemZone, BottomVegCover,season)) %>%
  select(which(!colSums(., na.rm=TRUE) %in% 0))

SXSf_winter_env <- SXS_full %>%
  filter(season == "winter") %>%
  subset(select = c(system,seasonYear,BottomVegCover))

winterfAbundPERM = adonis2(SXSf_winter_spe ~ system+seasonYear+BottomVegCover,
                          data = SXSf_winter_env,
                          method="bray",
                          permutations=999)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = SXSf_winter_spe ~ system + seasonYear + BottomVegCover, data = SXSf_winter_env, permutations = 999, method = "bray")
# Df SumOfSqs      R2       F Pr(>F)    
# system            3    256.4 0.08045 264.887  0.001 ***
# seasonYear       22     76.2 0.02389  10.727  0.001 ***
# BottomVegCover    1    107.6 0.03376 333.459  0.001 ***
# Residual       8513   2746.9 0.86189                   
# Total          8539   3187.1 1.00000                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# TBFull <- SiteXSpeciesFull %>%
#   filter(system == "TB") %>%
#   subset(select = -c(system)) %>%
#   column_to_rownames(var = "seasonYear")

SXSf_summer_spe <- SXS_full %>%
  filter(season == "summer") %>%
  subset(select = -c(season, system, seasonYear, Reference,systemZone, BottomVegCover,season)) %>%
  select(which(!colSums(., na.rm=TRUE) %in% 0))

SXSf_summer_env <- SXS_full %>%
  filter(season == "summer") %>%
  subset(select = c(system,seasonYear,BottomVegCover))

summerfAbundPERM = adonis2(SXSf_summer_spe ~ system+seasonYear+BottomVegCover,
                           data = SXSf_summer_env,
                           method="bray",
                           permutations=999)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = SXSf_summer_spe ~ system + seasonYear + BottomVegCover, data = SXSf_summer_env, permutations = 999, method = "bray")
# Df SumOfSqs      R2        F Pr(>F)    
# system            3   235.08 0.08217 288.9500  0.001 ***
# seasonYear       22    40.60 0.01419   6.8047  0.001 ***
# BottomVegCover    1   164.84 0.05762 607.8462  0.001 ***
# Residual       8925  2420.36 0.84602                    
# Total          8951  2860.88 1.00000                    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



systemLogic <- unique(SiteXSpeciesFull$system)

out <- list()

library(vegan)

for (i in 1:length(systemLogic)){
  tempDF <- SiteXSpeciesFull %>%
    filter(system == as.character(systemLogic[i])) %>%
    subset(select = -c(system)) %>%
    column_to_rownames(var = "seasonYear") 
  #select(which(!colSums(., na.rm=TRUE) %in% 0)) #technically unnecessary removal of totally absent taxa
  
  distmat <- as.matrix(vegdist(tempDF, method ="bray"))
  
  name <- as.character(systemLogic[i])
  
  out[[name]] <- distmat[,1]
}

adonis(distmat ~ system,
       data = SiteXSpeciesFull %>% filter(system=="TB"),
       permutations=99, 
       method = "bray")

yearLogic <- SiteXSpeciesFull %>%
  subset(select = c(system, seasonYear)) %>%
  unique()

lessgo <- stack(out) %>%
  rename(commsim = values) %>%
  rename(system = ind) %>%
  cbind(seasonYear = yearLogic$seasonYear) %>%
  #filter(commsim > 0) %>%
  mutate(season = "winter")



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