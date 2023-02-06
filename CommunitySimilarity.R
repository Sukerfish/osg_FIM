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
  subset(select = c(system,seasonYear,BottomVegCover)) %>%
  mutate(contYear = as.numeric(as.character(seasonYear)))

# bray curtis distance... square rooted - to make all eigenvalues positive per 
SXSf_winter_bray = vegdist(SXSf_winter_spe)
SXSf_winter_bray2 = SXSf_winter_bray^.5

homdisp1 = permutest(betadisper(SXSf_winter_bray2, SXSf_winter_env$system, type = "centroid"), 
                     pairwise = TRUE, permutations = 999, parallel = 6)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups       3  1.794 0.59793 157.91    999  0.001 ***
#   Residuals 8536 32.322 0.00379                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Pairwise comparisons:
#   (Observed p-value below diagonal, permuted p-value above diagonal)
# AP          CH          CK    TB
# AP              1.0000e-03  1.0000e-03 0.001
# CH  4.1915e-06              1.0000e-03 0.001
# CK  1.2255e-08  3.9532e-29             0.001
# TB  6.7991e-42 8.5214e-104  2.0450e-09      
homdisp2 = permutest(betadisper(SXSf_winter_bray2, SXSf_winter_env$seasonYear, type = "centroid"), 
                     pairwise = TRUE, permutations = 999, parallel = 6)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)    
# Groups      22  0.6312 0.0286900 10.218    999  0.001 ***
#   Residuals 8517 23.9133 0.0028077                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Pairwise comparisons:
#   (Observed p-value below diagonal, permuted p-value above diagonal)
#            2001       2002       2003       2004       2005       2006       2007       2008       2009       2010       2011       2012       2013       2014       2015       2016       2017       2018       2019       2020       1998       1999  2000
# 2001            1.0000e-03 2.0400e-01 7.7600e-01 1.0000e-03 1.0000e-03 7.0000e-03 1.0000e-03 3.0000e-03 5.2000e-02 9.3300e-01 1.0000e-03 1.0000e-03 1.2000e-02 2.0000e-03 1.5000e-02 1.0000e-03 1.5000e-02 1.0000e-03 1.0000e-03 6.8900e-01 8.3100e-01 0.452
# 2002 2.2542e-04            1.2000e-02 1.0000e-03 6.7500e-01 4.3300e-01 1.8700e-01 4.1600e-01 4.3700e-01 1.0000e-03 1.0000e-03 6.0600e-01 8.5800e-01 1.1700e-01 5.8200e-01 1.0100e-01 6.6300e-01 1.3000e-01 2.5400e-01 5.3800e-01 1.0000e-03 2.0000e-03 0.001
# 2003 2.0824e-01 7.1322e-03            2.8400e-01 1.1000e-02 1.0000e-03 9.0000e-02 1.0000e-03 4.8000e-02 2.0000e-03 1.5300e-01 2.0000e-03 6.0000e-03 1.7200e-01 2.9000e-02 2.4900e-01 3.0000e-03 1.8800e-01 1.0000e-03 1.0000e-03 1.2300e-01 1.8100e-01 0.052
# 2004 7.7651e-01 2.1936e-04 2.9000e-01            1.0000e-03 1.0000e-03 5.0000e-03 1.0000e-03 3.0000e-03 1.2000e-02 7.0400e-01 1.0000e-03 1.0000e-03 1.3000e-02 2.0000e-03 3.5000e-02 1.0000e-03 8.0000e-03 1.0000e-03 1.0000e-03 5.1700e-01 6.4400e-01 0.309
# 2005 2.3567e-04 6.7627e-01 1.1045e-02 2.5898e-04            1.8700e-01 3.3100e-01 1.7900e-01 6.3600e-01 1.0000e-03 1.0000e-03 3.3800e-01 5.1400e-01 2.2500e-01 8.4700e-01 1.6900e-01 3.4600e-01 1.9800e-01 1.0300e-01 2.7000e-01 1.0000e-03 1.0000e-03 0.001
# 2006 9.4493e-07 4.4374e-01 1.1646e-04 7.1250e-07 1.8521e-01            1.7000e-02 9.3900e-01 8.5000e-02 1.0000e-03 1.0000e-03 7.8900e-01 5.4400e-01 1.3000e-02 1.6600e-01 1.0000e-02 7.4700e-01 9.0000e-03 6.7100e-01 9.0100e-01 1.0000e-03 1.0000e-03 0.001
# 2007 3.3074e-03 1.8165e-01 9.1086e-02 4.4599e-03 3.1754e-01 1.6282e-02            2.6000e-02 6.1200e-01 1.0000e-03 4.0000e-03 4.9000e-02 1.1100e-01 7.5200e-01 4.6100e-01 6.2800e-01 6.5000e-02 7.1600e-01 1.1000e-02 3.1000e-02 1.0000e-03 6.0000e-03 0.001
# 2008 2.1799e-06 4.3305e-01 2.0101e-04 1.5550e-06 1.8779e-01 9.4216e-01 1.9458e-02            1.0000e-01 1.0000e-03 1.0000e-03 7.4300e-01 5.3500e-01 1.1000e-02 1.7200e-01 1.0000e-02 6.8800e-01 1.6000e-02 7.4800e-01 8.6300e-01 1.0000e-03 1.0000e-03 0.001
# 2009 1.9426e-03 4.3243e-01 4.7830e-02 2.2683e-03 6.6235e-01 9.3776e-02 6.2868e-01 9.7281e-02            1.0000e-03 1.0000e-03 1.6600e-01 3.2300e-01 4.6200e-01 8.0900e-01 3.8800e-01 2.1700e-01 4.2900e-01 4.3000e-02 1.5200e-01 1.0000e-03 3.0000e-03 0.001
# 2010 4.5083e-02 3.8191e-09 5.2988e-04 1.2909e-02 1.1247e-09 2.1818e-13 6.7274e-08 8.1735e-13 5.1293e-08            3.5000e-02 1.0000e-03 1.0000e-03 1.0000e-03 1.0000e-03 1.0000e-03 1.0000e-03 1.0000e-03 1.0000e-03 1.0000e-03 1.5400e-01 1.4300e-01 0.323
# 2011 9.2763e-01 6.7725e-05 1.5200e-01 6.8743e-01 6.8142e-05 1.6183e-07 1.3347e-03 3.4077e-07 6.8358e-04 3.8455e-02            1.0000e-03 1.0000e-03 8.0000e-03 1.0000e-03 1.0000e-02 1.0000e-03 6.0000e-03 1.0000e-03 1.0000e-03 7.5700e-01 8.9600e-01 0.491
# 2012 1.0148e-05 6.2305e-01 7.2501e-04 8.1128e-06 3.2325e-01 8.0206e-01 4.7302e-02 7.6115e-01 1.7786e-01 9.1699e-12 1.8640e-06            7.5400e-01 2.9000e-02 2.6200e-01 2.2000e-02 9.4300e-01 2.8000e-02 5.3000e-01 9.0800e-01 1.0000e-03 1.0000e-03 0.001
# 2013 4.8060e-05 8.4643e-01 2.5630e-03 4.2229e-05 5.1666e-01 5.6851e-01 1.0694e-01 5.4408e-01 3.0540e-01 1.0191e-10 1.0061e-05 7.5975e-01            6.3000e-02 4.1100e-01 5.3000e-02 8.1200e-01 6.6000e-02 3.4300e-01 6.7500e-01 1.0000e-03 1.0000e-03 0.001
# 2014 9.9319e-03 1.2145e-01 1.7886e-01 1.3489e-02 2.0926e-01 9.5223e-03 7.5755e-01 1.1409e-02 4.5463e-01 5.7984e-07 4.4684e-03 2.8428e-02 6.6487e-02            3.1100e-01 8.7100e-01 3.6000e-02 9.8300e-01 5.0000e-03 2.4000e-02 3.0000e-03 1.3000e-02 0.003
# 2015 8.1440e-04 5.7223e-01 2.5502e-02 9.0099e-04 8.4981e-01 1.5217e-01 4.5472e-01 1.5348e-01 8.1257e-01 9.7785e-09 2.4926e-04 2.6411e-01 4.2609e-01 3.1488e-01            2.6300e-01 2.7800e-01 3.1800e-01 7.1000e-02 2.3400e-01 2.0000e-03 2.0000e-03 0.001
# 2016 1.7691e-02 1.0091e-01 2.4637e-01 2.3965e-02 1.7027e-01 7.6983e-03 6.4619e-01 9.2360e-03 3.8285e-01 2.1558e-06 8.5954e-03 2.2929e-02 5.3725e-02 8.7740e-01 2.6095e-01            2.8000e-02 8.9600e-01 3.0000e-03 1.6000e-02 1.3000e-02 2.5000e-02 0.005
# 2017 1.6377e-05 6.7662e-01 1.0767e-03 1.3203e-05 3.6872e-01 7.4739e-01 6.0357e-02 7.1011e-01 2.0641e-01 1.6325e-11 2.9443e-06 9.4307e-01 8.1637e-01 3.6350e-02 3.0075e-01 2.9228e-02            4.0000e-02 4.9500e-01 8.3200e-01 1.0000e-03 1.0000e-03 0.001
# 2018 1.0614e-02 1.1678e-01 1.8622e-01 1.4458e-02 2.0139e-01 8.8910e-03 7.3999e-01 1.0750e-02 4.4240e-01 6.8459e-07 4.8495e-03 2.6986e-02 6.3632e-02 9.8208e-01 3.0520e-01 8.9477e-01 3.4637e-02            4.0000e-03 2.7000e-02 6.0000e-03 2.1000e-02 0.003
# 2019 2.5578e-07 2.6714e-01 3.4178e-05 1.6144e-07 9.0450e-02 6.7599e-01 5.8137e-03 7.4849e-01 4.4019e-02 3.0146e-14 3.4419e-08 5.2849e-01 3.5043e-01 3.3648e-03 7.5333e-02 2.7593e-03 4.8776e-01 3.1310e-03            5.7900e-01 1.0000e-03 1.0000e-03 0.001
# 2020 5.3257e-06 5.4683e-01 4.3229e-04 4.0812e-06 2.6503e-01 8.9539e-01 3.3727e-02 8.4776e-01 1.4191e-01 3.1476e-12 9.1205e-07 9.0981e-01 6.7505e-01 1.9996e-02 2.1606e-01 1.6117e-02 8.5427e-01 1.8923e-02 6.0468e-01            1.0000e-03 1.0000e-03 0.001
# 1998 7.0941e-01 1.3827e-04 1.2049e-01 5.1249e-01 1.2343e-04 3.4646e-07 1.5687e-03 1.5520e-06 1.3589e-03 1.7666e-01 7.5682e-01 7.5525e-06 4.1270e-05 5.6703e-03 6.0507e-04 1.0765e-02 1.4960e-05 5.9774e-03 1.3469e-07 3.8690e-06            8.8600e-01 0.728
# 1999 8.2817e-01 4.2137e-04 1.7795e-01 6.3038e-01 4.0287e-04 2.4037e-06 3.9744e-03 7.5933e-06 3.1229e-03 1.3921e-01 8.8078e-01 3.0587e-05 1.3315e-04 1.1786e-02 1.5009e-03 2.0421e-02 5.3286e-05 1.2395e-02 9.3907e-07 1.7007e-05 8.9144e-01            0.663
# 2000 4.8051e-01 5.4754e-05 6.0538e-02 3.1277e-01 4.1395e-05 1.2306e-07 5.3748e-04 5.0484e-07 4.9104e-04 3.2971e-01 5.0654e-01 2.4578e-06 1.3342e-05 2.0829e-03 2.0741e-04 4.2121e-03 4.6789e-06 2.2146e-03 4.6001e-08 1.2489e-06 7.4956e-01 6.6108e-01          

winterfAbundPERM = adonis2(SXSf_winter_bray2 ~ system * seasonYear + contYear + BottomVegCover,
                          data = SXSf_winter_env,
                          #add = "lingoes",
                          parallel = 6,
                          #method="bray", 
                          permutations=999)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = SXSf_winter_bray2 ~ system * seasonYear + contYear + BottomVegCover, data = SXSf_winter_env, permutations = 999, parallel = 6)
#                     Df SumOfSqs      R2        F Pr(>F)    
# system               3    170.0 0.04687 146.4418  0.001 ***
# seasonYear          22     55.8 0.01539   6.5585  0.001 ***
# BottomVegCover       1     71.2 0.01963 183.9698  0.001 ***
# system:seasonYear   60     59.0 0.01628   2.5430  0.001 ***
# Residual          8453   3270.9 0.90183                    
# Total             8539   3627.0 1.00000                    
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
  subset(select = c(system,seasonYear,BottomVegCover)) %>%
  mutate(contYear = as.numeric(as.character(seasonYear)))

# bray curtis distance... square rooted - to make all eigenvalues positive per 
SXSf_summer_bray = vegdist(SXSf_summer_spe)
SXSf_summer_bray2 = SXSf_summer_bray^.5

homdispSySummer = permutest(betadisper(SXSf_summer_bray2, SXSf_summer_env$system, type = "centroid"), 
                     pairwise = TRUE, permutations = 999, parallel = 6)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#             Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups       3  3.6404 1.21346 348.41    999  0.001 ***
# Residuals 8948 31.1645 0.00348                         
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Pairwise comparisons:
# (Observed p-value below diagonal, permuted p-value above diagonal)
# AP          CH          CK    TB
# AP              1.0000e-03  5.0000e-03 0.054
# CH 2.9358e-110              1.0000e-03 0.001
# CK  1.8808e-03 3.0944e-141             0.001
# TB  4.3747e-02 6.2034e-121  1.2539e-07  

homdispYrSummer = permutest(betadisper(SXSf_winter_bray2, SXSf_winter_env$seasonYear, type = "centroid"), 
                     pairwise = TRUE, permutations = 999, parallel = 6)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#             Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)    
# Groups      22  0.6312 0.0286900 10.218    999  0.001 ***
# Residuals 8517 23.9133 0.0028077                         
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Pairwise comparisons:
# (Observed p-value below diagonal, permuted p-value above diagonal)
# 2001       2002       2003       2004       2005       2006       2007       2008       2009       2010       2011       2012       2013       2014       2015       2016
# 2001            1.0000e-03 1.9600e-01 7.7700e-01 1.0000e-03 1.0000e-03 2.0000e-03 1.0000e-03 3.0000e-03 4.2000e-02 9.3700e-01 1.0000e-03 1.0000e-03 7.0000e-03 1.0000e-03 1.7000e-02
# 2002 2.2542e-04            6.0000e-03 1.0000e-03 6.7700e-01 4.3200e-01 1.8900e-01 4.1300e-01 4.5400e-01 1.0000e-03 1.0000e-03 5.9900e-01 8.4300e-01 1.1100e-01 5.5000e-01 8.8000e-02
# 2003 2.0824e-01 7.1322e-03            2.4100e-01 9.0000e-03 2.0000e-03 8.8000e-02 2.0000e-03 4.9000e-02 2.0000e-03 1.3300e-01 1.0000e-03 4.0000e-03 1.8000e-01 2.8000e-02 2.5700e-01
# 2004 7.7651e-01 2.1936e-04 2.9000e-01            1.0000e-03 1.0000e-03 6.0000e-03 1.0000e-03 2.0000e-03 1.1000e-02 6.7700e-01 1.0000e-03 1.0000e-03 1.5000e-02 1.0000e-03 2.4000e-02
# 2005 2.3567e-04 6.7627e-01 1.1045e-02 2.5898e-04            1.9400e-01 3.3300e-01 1.8400e-01 6.8400e-01 1.0000e-03 1.0000e-03 3.1600e-01 5.3900e-01 2.0900e-01 8.4200e-01 1.6000e-01
# 2006 9.4493e-07 4.4374e-01 1.1646e-04 7.1250e-07 1.8521e-01            2.0000e-02 9.4200e-01 9.4000e-02 1.0000e-03 1.0000e-03 8.1800e-01 5.8200e-01 8.0000e-03 1.6100e-01 8.0000e-03
# 2007 3.3074e-03 1.8165e-01 9.1086e-02 4.4599e-03 3.1754e-01 1.6282e-02            2.2000e-02 6.4100e-01 1.0000e-03 3.0000e-03 4.4000e-02 1.1200e-01 7.7900e-01 4.7500e-01 6.1900e-01
# 2008 2.1799e-06 4.3305e-01 2.0101e-04 1.5550e-06 1.8779e-01 9.4216e-01 1.9458e-02            9.1000e-02 1.0000e-03 1.0000e-03 7.5500e-01 5.5400e-01 1.4000e-02 1.4600e-01 9.0000e-03
# 2009 1.9426e-03 4.3243e-01 4.7830e-02 2.2683e-03 6.6235e-01 9.3776e-02 6.2868e-01 9.7281e-02            1.0000e-03 1.0000e-03 1.8700e-01 3.0900e-01 4.4600e-01 8.2700e-01 3.5900e-01
# 2010 4.5083e-02 3.8191e-09 5.2988e-04 1.2909e-02 1.1247e-09 2.1818e-13 6.7274e-08 8.1735e-13 5.1293e-08            4.3000e-02 1.0000e-03 1.0000e-03 1.0000e-03 1.0000e-03 1.0000e-03
# 2011 9.2763e-01 6.7725e-05 1.5200e-01 6.8743e-01 6.8142e-05 1.6183e-07 1.3347e-03 3.4077e-07 6.8358e-04 3.8455e-02            1.0000e-03 1.0000e-03 5.0000e-03 1.0000e-03 9.0000e-03
# 2012 1.0148e-05 6.2305e-01 7.2501e-04 8.1128e-06 3.2325e-01 8.0206e-01 4.7302e-02 7.6115e-01 1.7786e-01 9.1699e-12 1.8640e-06            7.4400e-01 2.9000e-02 2.5100e-01 2.1000e-02
# 2013 4.8060e-05 8.4643e-01 2.5630e-03 4.2229e-05 5.1666e-01 5.6851e-01 1.0694e-01 5.4408e-01 3.0540e-01 1.0191e-10 1.0061e-05 7.5975e-01            6.8000e-02 4.4200e-01 5.1000e-02
# 2014 9.9319e-03 1.2145e-01 1.7886e-01 1.3489e-02 2.0926e-01 9.5223e-03 7.5755e-01 1.1409e-02 4.5463e-01 5.7984e-07 4.4684e-03 2.8428e-02 6.6487e-02            2.9900e-01 8.6900e-01
# 2015 8.1440e-04 5.7223e-01 2.5502e-02 9.0099e-04 8.4981e-01 1.5217e-01 4.5472e-01 1.5348e-01 8.1257e-01 9.7785e-09 2.4926e-04 2.6411e-01 4.2609e-01 3.1488e-01            2.4300e-01
# 2016 1.7691e-02 1.0091e-01 2.4637e-01 2.3965e-02 1.7027e-01 7.6983e-03 6.4619e-01 9.2360e-03 3.8285e-01 2.1558e-06 8.5954e-03 2.2929e-02 5.3725e-02 8.7740e-01 2.6095e-01           
# 2017 1.6377e-05 6.7662e-01 1.0767e-03 1.3203e-05 3.6872e-01 7.4739e-01 6.0357e-02 7.1011e-01 2.0641e-01 1.6325e-11 2.9443e-06 9.4307e-01 8.1637e-01 3.6350e-02 3.0075e-01 2.9228e-02
# 2018 1.0614e-02 1.1678e-01 1.8622e-01 1.4458e-02 2.0139e-01 8.8910e-03 7.3999e-01 1.0750e-02 4.4240e-01 6.8459e-07 4.8495e-03 2.6986e-02 6.3632e-02 9.8208e-01 3.0520e-01 8.9477e-01
# 2019 2.5578e-07 2.6714e-01 3.4178e-05 1.6144e-07 9.0450e-02 6.7599e-01 5.8137e-03 7.4849e-01 4.4019e-02 3.0146e-14 3.4419e-08 5.2849e-01 3.5043e-01 3.3648e-03 7.5333e-02 2.7593e-03
# 2020 5.3257e-06 5.4683e-01 4.3229e-04 4.0812e-06 2.6503e-01 8.9539e-01 3.3727e-02 8.4776e-01 1.4191e-01 3.1476e-12 9.1205e-07 9.0981e-01 6.7505e-01 1.9996e-02 2.1606e-01 1.6117e-02
# 1998 7.0941e-01 1.3827e-04 1.2049e-01 5.1249e-01 1.2343e-04 3.4646e-07 1.5687e-03 1.5520e-06 1.3589e-03 1.7666e-01 7.5682e-01 7.5525e-06 4.1270e-05 5.6703e-03 6.0507e-04 1.0765e-02
# 1999 8.2817e-01 4.2137e-04 1.7795e-01 6.3038e-01 4.0287e-04 2.4037e-06 3.9744e-03 7.5933e-06 3.1229e-03 1.3921e-01 8.8078e-01 3.0587e-05 1.3315e-04 1.1786e-02 1.5009e-03 2.0421e-02
# 2000 4.8051e-01 5.4754e-05 6.0538e-02 3.1277e-01 4.1395e-05 1.2306e-07 5.3748e-04 5.0484e-07 4.9104e-04 3.2971e-01 5.0654e-01 2.4578e-06 1.3342e-05 2.0829e-03 2.0741e-04 4.2121e-03
# 2017       2018       2019       2020       1998       1999  2000
# 2001 1.0000e-03 1.5000e-02 1.0000e-03 1.0000e-03 7.2600e-01 8.1500e-01 0.463
# 2002 6.8300e-01 1.2400e-01 2.6600e-01 5.8700e-01 2.0000e-03 2.0000e-03 0.001
# 2003 1.0000e-03 1.6300e-01 1.0000e-03 3.0000e-03 1.2200e-01 1.7900e-01 0.057
# 2004 1.0000e-03 1.9000e-02 1.0000e-03 1.0000e-03 5.2400e-01 6.3100e-01 0.318
# 2005 3.6600e-01 2.1300e-01 9.4000e-02 2.8800e-01 1.0000e-03 1.0000e-03 0.001
# 2006 7.2700e-01 8.0000e-03 6.7300e-01 8.8900e-01 1.0000e-03 1.0000e-03 0.001
# 2007 6.7000e-02 7.4400e-01 6.0000e-03 4.4000e-02 4.0000e-03 4.0000e-03 0.001
# 2008 7.0900e-01 1.5000e-02 7.3600e-01 8.5800e-01 1.0000e-03 1.0000e-03 0.001
# 2009 2.0000e-01 4.2800e-01 3.4000e-02 1.3900e-01 2.0000e-03 5.0000e-03 0.002
# 2010 1.0000e-03 1.0000e-03 1.0000e-03 1.0000e-03 1.7600e-01 1.5600e-01 0.317
# 2011 1.0000e-03 4.0000e-03 1.0000e-03 1.0000e-03 7.3500e-01 9.0400e-01 0.529
# 2012 9.3200e-01 2.8000e-02 5.1900e-01 9.2200e-01 1.0000e-03 1.0000e-03 0.001
# 2013 8.2000e-01 6.7000e-02 3.4700e-01 6.8000e-01 1.0000e-03 2.0000e-03 0.001
# 2014 4.1000e-02 9.7400e-01 4.0000e-03 1.8000e-02 4.0000e-03 1.8000e-02 0.003
# 2015 3.0200e-01 3.0300e-01 7.3000e-02 2.1100e-01 2.0000e-03 3.0000e-03 0.001
# 2016 2.9000e-02 8.9600e-01 2.0000e-03 2.5000e-02 1.4000e-02 3.0000e-02 0.006
# 2017            3.3000e-02 4.4700e-01 8.5600e-01 1.0000e-03 1.0000e-03 0.001
# 2018 3.4637e-02            1.0000e-03 1.7000e-02 8.0000e-03 1.9000e-02 0.003
# 2019 4.8776e-01 3.1310e-03            5.9800e-01 1.0000e-03 1.0000e-03 0.001
# 2020 8.5427e-01 1.8923e-02 6.0468e-01            1.0000e-03 1.0000e-03 0.001
# 1998 1.4960e-05 5.9774e-03 1.3469e-07 3.8690e-06            9.1400e-01 0.742
# 1999 5.3286e-05 1.2395e-02 9.3907e-07 1.7007e-05 8.9144e-01            0.643
# 2000 4.6789e-06 2.2146e-03 4.6001e-08 1.2489e-06 7.4956e-01 6.6108e-01      

summerfAbundPERM = adonis2(SXSf_summer_spe ~ system * seasonYear + contYear + BottomVegCover,
                           data = SXSf_summer_env,
                           #add = "lingoes",
                           parallel = 6,
                           #method="bray", 
                           permutations=999)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = SXSf_summer_spe ~ system * seasonYear + contYear + BottomVegCover, data = SXSf_summer_env, permutations = 999, parallel = 6)
#                     Df SumOfSqs      R2        F Pr(>F)    
# system               3   235.08 0.08217 293.0915  0.001 ***
# seasonYear          22    40.60 0.01419   6.9022  0.001 ***
# BottomVegCover       1   164.84 0.05762 616.5584  0.001 ***
# system:seasonYear   60    50.24 0.01756   3.1320  0.001 ***
# Residual          8865  2370.12 0.82846                    
# Total             8951  2860.88 1.00000                    
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