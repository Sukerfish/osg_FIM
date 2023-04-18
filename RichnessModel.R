# #### Import and Filter Data #####

library(tidyverse)

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

SXS_filteredList <- list()
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
    mutate(n = sum(across(!c(Reference:seasonYear)))) %>%
    ungroup() %>%
    filter(n > 0) %>%
    select(Reference)
 
SXS_filteredList[[i]] <- df_filtered
  
}

SXS_filtered <- bind_rows(SXS_filteredList) %>%
  left_join(SXS_full)

SXS_filtered_env <- SXS_filtered %>%
  select(c(Reference, systemSeason, seasonYear, systemZone, BottomVegCover, Temperature))

SXS_filtered_spp <- SXS_filtered %>%
  select(!c(Reference, systemSeason, seasonYear, systemZone, BottomVegCover, Temperature))

SXR_filtered_spp <- SXS_filtered %>%
  select(!c(systemSeason, seasonYear, systemZone, BottomVegCover, Temperature)) %>%
  pivot_longer(cols = !c(Reference),
               names_to = "taxa") %>%
  mutate(PA = ifelse(value > 0, 1, 0)) %>%
  select(!c(taxa, value)) %>%
  group_by(Reference) %>%
  mutate(N = sum(PA)) %>%
  select(!c(PA)) %>%
  distinct()

SXR_filtered <- SXR_filtered_spp %>%
  left_join(SXS_filtered_env) %>%
  group_by(systemSeason) %>%
  add_count(name = "n_hauls") %>%
  ungroup()

modelDFM <- SXR_filtered %>%
  separate(systemSeason,
           c("system","season"),
           sep = "_") %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH"))) %>%
  mutate(contYear = as.numeric(as.character(seasonYear))) %>%
  mutate(seasonYear = as.factor(seasonYear)) %>%
  # mutate(n_hauls = log(n_hauls)) %>%
  # filter(season == "summer") %>%
  # filter(system == "TB") %>%
  filter(season == "winter")

modelDF <- as.data.frame(modelDFM)

# modelDF_summer <- FullRichness %>%
#   group_by(seasonYear, system, season) %>%
#   add_count(name = "n_hauls") %>%
#   ungroup() %>%
#   mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH"))) %>%
#   mutate(contYear = as.numeric(as.character(seasonYear))) %>%
#   mutate(seasonYear = as.factor(seasonYear)) %>%
#   # mutate(n_hauls = log(n_hauls)) %>%
#   # filter(season == "summer") %>%
#   # filter(system == "TB") %>%
#   filter(season == "summer")


# bottomVeg <- TidyBio %>%
#   subset(select = c(Reference, BottomVegCover)) %>%
#   group_by(Reference) %>%
#   summarise(bvc.mean = mean(BottomVegCover, na.rm = TRUE)) %>%
#   mutate(BottomVegCover = bvc.mean) %>%
#   filter(!is.na(bvc.mean)) %>%
#   subset(select = c(Reference, BottomVegCover))
# 
# OnlyFish <- TidyBio %>%
#   filter(Scientificname != "No fish")
# 
# NoFish <- TidyBio %>%
#   filter(Scientificname == "No fish") %>%
#   mutate(n = 0) %>%
#   subset(select = c(Reference, n, system, season, seasonYear, systemZone))
# 
# FullRichness <- OnlyFish %>%
#   group_by(Reference) %>%
#   count() %>%
#   inner_join(TidyRefsList) %>%
#   ungroup() %>%
#   bind_rows(NoFish) %>%
#   left_join(bottomVeg) 



# GLMM:
# Response: richness or total abundance
# Fixed effects: estuary (four levels) and year
# Random effect: within-estuary zone
# Offset term for variation in number of seine hauls
# Autoregressive term to account for serial correlation



ggplot(modelDF_winter, aes(x=n)) +
  geom_histogram(binwidth=1) +
  theme(axis.text=element_text(size = 12)) +
  theme(axis.title=element_text(size = 16)) +
  theme(strip.text = element_text(size = 16)) +
  ggtitle("Richness per haul") +
  theme(title=element_text(size = 20)) +
  facet_grid(season~system)

ggplot(modelDF_winter, aes(y=n,
                    x = factor(seasonYear))) +
  geom_boxplot() +
  scale_x_discrete(breaks=c(2000,2005,2010,2015,2020))+
  theme(axis.text=element_text(size = 12)) +
  theme(axis.title=element_text(size = 16)) +
  theme(strip.text = element_text(size = 16)) +
  ggtitle("Richness per haul over time") +
  theme(title=element_text(size = 20)) +
  facet_grid(season~system)


# test <-modelDF %>%
#   subset(select = c(seasonYear))
# df <- fct_reorder(factor(modelDF$seasonYear),test$seasonYear)

###### pql ####
library(nlme)
library(lme4)
library(MASS) # needs MASS (version 7.3-58)

#DOES NOT WORK
glmmPQL <- glmmPQL(N ~ system + seasonYear + offset(log(n_hauls)),
                    random = ~ 1|systemZone,
                    family = poisson,
                    #correlation = corARMA(form = ~1|systemZone, p = 1, q = 1),
                    data = modelDF)
summary(glmmPQL)
plot(glmmPQL)

#alternative cor structure ... corAR1 converges
glmmPQL <- glmmPQL(n ~ system + seasonYear + offset(log(n_hauls)),
                   random = ~ 1|systemZone,
                   family = poisson,
                   correlation = corAR1(form = ~1|systemZone),
                   data = modelDF)
summary(glmmPQL)
plot(glmmPQL)

#converges without correlation
glmmPQL <- glmmPQL(n ~ system + seasonYear + offset(log(n_hauls)),
                   random = ~ 1|systemZone,
                   family = poisson,
                   data = modelDF)
summary(glmmPQL)
plot(glmmPQL)

testlmer <- glmer(n ~ system + seasonYear + offset(log(n_hauls))
                   + (1|systemZone),
                   family = poisson,
                   #correlation = corAR1(form = ~1|systemZone),
                   data = modelDF)
summary(testlmer)



##### glmmTMB ####
library(glmmTMB)
rmod_tmb <- glmmTMB(N~system+
                      seasonYear+
                      offset(log(n_hauls))+
                      #ar1(factor(seasonYear))+
                      (1|systemZone)+
                      (1|Temperature)+
                      (1|BottomVegCover),
                    zi=~0,
                    family=poisson,
                    data=modelDF)
summary(rmod_tmb)
VarCorr(rmod_tmb)

##### lme4 ####
#converges if seasonYear is random
library(lme4)
library(lmerTest)

glme_rich <- glmer(N~system+
                       seasonYear+
                       offset(log(n_hauls))+
                       (1|BottomVegCover)+
                       (1|systemZone)+
                     (1|Temperature)+
                       (1|contYear),
                     family="poisson",
                     data=modelDF,
                     control=glmerControl(optimizer="bobyqa",
                                          check.conv.grad=.makeCC("warning",0.05)))

summary(glme_rich)
#ranef(glme_rich)

# WHEN contYear is included:
# boundary (singular) fit: see help('isSingular')
# 
# > ranef(gmod_lme4_L)
# $BottomVegCover
# (Intercept)
# 0   -0.286778197
# 1   -0.554656884
# 2   -0.485282983
# 3   -0.316112389
# 4   -0.153897614
# 5   -0.420456838
# 6    0.052852022
# 7   -0.095967185
# 8   -0.041994018
# 9    0.182975981
# 10  -0.277217398
# 12  -0.333955075
# 13   0.199986612
# 15  -0.141748694
# 17  -0.225760348
# 20  -0.094252424
# 24  -0.338278760
# 25  -0.019703572
# 28   0.063680466
# 30   0.038523397
# 35   0.044137427
# 38  -0.004546357
# 40   0.056430271
# 45   0.142677940
# 50   0.139696794
# 55   0.333678760
# 60   0.201555754
# 65   0.118314700
# 69   0.248448168
# 70   0.217038225
# 73   0.039055369
# 75   0.174176325
# 78   0.236469920
# 79   0.037584506
# 80   0.214267889
# 85   0.291731016
# 90   0.229240836
# 92   0.002031692
# 94   0.050437252
# 95   0.198828684
# 96   0.071597185
# 97  -0.213756871
# 98   0.117693689
# 99   0.381708174
# 100  0.089781643
# 101 -0.005640213
# 
# $contYear
# (Intercept)
# 1998           0
# 1999           0
# 2000           0
# 2001           0
# 2002           0
# 2003           0
# 2004           0
# 2005           0
# 2006           0
# 2007           0
# 2008           0
# 2009           0
# 2010           0
# 2011           0
# 2012           0
# 2013           0
# 2014           0
# 2015           0
# 2016           0
# 2017           0
# 2018           0
# 2019           0
# 2020           0
# 
# $systemZone
# (Intercept)
# AP_A  0.125793430
# AP_B -0.123653260
# CH_A  0.056074224
# CH_B  0.011136499
# CH_C  0.047730151
# CH_D -0.110828474
# CK_B -0.021993325
# CK_C  0.024001857
# TB_A  0.105474887
# TB_B -0.046134171
# TB_C -0.019111862
# TB_D -0.043769500
# TB_E  0.008875772
# 
# with conditional variances for “BottomVegCover” “contYear” “systemZone” 
summary(glme_rich_summer)
plot(glme_rich_summer)
plot(glme_rich_summer,systemZone~resid(.))


glme_rich_summer <- glmer(n~system+
                       seasonYear+
                       offset(log(n_hauls))+
                       (1|BottomVegCover)+
                       (1|systemZone),
                     family="poisson",
                     data=modelDF_summer,
                     control=glmerControl(optimizer="bobyqa",
                                          check.conv.grad=.makeCC("warning",0.05)))
ranef(glme_rich_summer)
summary(glme_rich_summer)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: n ~ system + seasonYear + offset(log(n_hauls)) + (1 | BottomVegCover) +  
#   (1 | systemZone)
# Data: modelDF_summer
# Control: glmerControl(optimizer = "bobyqa", check.conv.grad = .makeCC("warning",      0.05))
# 
# AIC      BIC   logLik deviance df.resid 
# 53437.9  53636.9 -26691.0  53381.9     8979 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -3.2559 -1.0306 -0.1251  0.8969  6.7312 
# 
# Random effects:
#   Groups         Name        Variance Std.Dev.
# BottomVegCover (Intercept) 0.067345 0.25951 
# systemZone     (Intercept) 0.005563 0.07458 
# Number of obs: 9007, groups:  BottomVegCover, 46; systemZone, 13
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -1.714599   0.073560 -23.309  < 2e-16 ***
#   systemCK       -0.072798   0.075572  -0.963    0.335    
# systemTB       -0.660035   0.063317 -10.424  < 2e-16 ***
#   systemCH       -0.666871   0.065484 -10.184  < 2e-16 ***
#   seasonYear1999  0.010660   0.033726   0.316    0.752    
# seasonYear2000 -0.002186   0.033806  -0.065    0.948    
# seasonYear2001 -0.168063   0.029977  -5.606 2.07e-08 ***
#   seasonYear2002 -0.178254   0.030073  -5.927 3.08e-09 ***
#   seasonYear2003 -0.138177   0.029717  -4.650 3.32e-06 ***
#   seasonYear2004 -0.215828   0.029101  -7.416 1.20e-13 ***
#   seasonYear2005 -0.424016   0.029437 -14.404  < 2e-16 ***
#   seasonYear2006 -0.378166   0.029125 -12.984  < 2e-16 ***
#   seasonYear2007 -0.392300   0.029232 -13.420  < 2e-16 ***
#   seasonYear2008 -0.338052   0.028908 -11.694  < 2e-16 ***
#   seasonYear2009 -0.366478   0.029043 -12.618  < 2e-16 ***
#   seasonYear2010 -0.371903   0.029026 -12.813  < 2e-16 ***
#   seasonYear2011 -0.373310   0.028999 -12.873  < 2e-16 ***
#   seasonYear2012 -0.438333   0.029357 -14.931  < 2e-16 ***
#   seasonYear2013 -0.409770   0.029196 -14.035  < 2e-16 ***
#   seasonYear2014 -0.368040   0.028956 -12.711  < 2e-16 ***
#   seasonYear2015 -0.396270   0.029029 -13.651  < 2e-16 ***
#   seasonYear2016 -0.343117   0.028871 -11.884  < 2e-16 ***
#   seasonYear2017 -0.350129   0.028993 -12.076  < 2e-16 ***
#   seasonYear2018 -0.379346   0.029058 -13.055  < 2e-16 ***
#   seasonYear2019 -0.341449   0.028863 -11.830  < 2e-16 ***
#   seasonYear2020 -0.355158   0.028927 -12.278  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

glme_rich_winter <- glmer(n~system+
                            seasonYear+
                            offset(log(n_hauls))+
                            (1|BottomVegCover)+
                            (1|systemZone),
                          family="poisson",
                          data=modelDF_winter,
                          control=glmerControl(optimizer="bobyqa",
                                               check.conv.grad=.makeCC("warning",0.05)))
ranef(glme_rich_winter)
summary(glme_rich_winter)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: n ~ system + seasonYear + offset(log(n_hauls)) + (1 | BottomVegCover) +  
#   (1 | systemZone)
# Data: modelDF_winter
# Control: glmerControl(optimizer = "bobyqa", check.conv.grad = .makeCC("warning",      0.05))
# 
# AIC      BIC   logLik deviance df.resid 
# 47522.8  47721.7 -23733.4  47466.8     8950 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.9616 -1.0542 -0.1665  0.8835  6.9494 
# 
# Random effects:
#   Groups         Name        Variance Std.Dev.
# BottomVegCover (Intercept) 0.04696  0.2167  
# systemZone     (Intercept) 0.01988  0.1410  
# Number of obs: 8978, groups:  BottomVegCover, 47; systemZone, 13
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -2.524486   0.112436 -22.453  < 2e-16 ***
#   systemCK       -0.262329   0.142140  -1.846 0.064954 .  
# systemTB       -0.424898   0.118906  -3.573 0.000352 ***
#   systemCH       -0.158834   0.122977  -1.292 0.196503    
# seasonYear1999  0.004921   0.041893   0.117 0.906493    
# seasonYear2000  0.129259   0.040912   3.159 0.001581 ** 
#   seasonYear2001  0.067912   0.038215   1.777 0.075554 .  
# seasonYear2002 -0.155348   0.038880  -3.996 6.45e-05 ***
#   seasonYear2003  0.060321   0.037483   1.609 0.107549    
# seasonYear2004 -0.117618   0.036794  -3.197 0.001390 ** 
#   seasonYear2005 -0.333318   0.037394  -8.914  < 2e-16 ***
#   seasonYear2006 -0.344907   0.037255  -9.258  < 2e-16 ***
#   seasonYear2007 -0.229335   0.036449  -6.292 3.14e-10 ***
#   seasonYear2008 -0.351031   0.037212  -9.433  < 2e-16 ***
#   seasonYear2009 -0.454898   0.037774 -12.043  < 2e-16 ***
#   seasonYear2010 -0.259740   0.036674  -7.082 1.42e-12 ***
#   seasonYear2011 -0.333853   0.037093  -9.000  < 2e-16 ***
#   seasonYear2012 -0.398587   0.037512 -10.626  < 2e-16 ***
#   seasonYear2013 -0.277981   0.036759  -7.562 3.96e-14 ***
#   seasonYear2014 -0.198911   0.036336  -5.474 4.39e-08 ***
#   seasonYear2015 -0.269293   0.036613  -7.355 1.91e-13 ***
#   seasonYear2016 -0.224745   0.036466  -6.163 7.13e-10 ***
#   seasonYear2017 -0.263274   0.036618  -7.190 6.49e-13 ***
#   seasonYear2018 -0.231000   0.036502  -6.328 2.48e-10 ***
#   seasonYear2019 -0.283648   0.036705  -7.728 1.10e-14 ***
#   seasonYear2020 -0.294739   0.036857  -7.997 1.28e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


interceptOnly <- lmer(n ~ 1 + (1|systemZone),
                       data = modelDF)
summary(interceptOnly)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: n ~ 1 + (1 | systemZone)
#    Data: modelDF
# 
# REML criterion at convergence: 47233.5
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -2.3774 -0.7433 -0.1415  0.6357  4.8244 
# 
# Random effects:
#  Groups     Name        Variance Std.Dev.
#  systemZone (Intercept)  2.498   1.581   
#  Residual               11.045   3.323   
# Number of obs: 9002, groups:  systemZone, 13
# 
# Fixed effects:
#             Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)   5.5275     0.4398 12.0087   12.57 2.86e-08 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

model1 <- lmer(n ~ 1 + system + seasonYear + (1|systemZone),
               data = modelDF)
summary(model1)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: n ~ 1 + system + seasonYear + (1 | systemZone)
#    Data: modelDF
# 
# REML criterion at convergence: 47224.1
# 
# Scaled residuals:
#     Min      1Q  Median      3Q     Max
# -2.3921 -0.7530 -0.1404  0.6483  4.8387
# 
# Random effects:
#  Groups     Name        Variance Std.Dev.
#  systemZone (Intercept)  1.045   1.022
#  Residual               11.045   3.323
# Number of obs: 9002, groups:  systemZone, 13
# 
# Fixed effects:
#               Estimate Std. Error         df t value Pr(>|t|)
# (Intercept) -4.350e+00  1.149e+01  8.906e+03  -0.379    0.705
# systemCK    -1.297e+00  1.029e+00  8.945e+00  -1.261    0.239
# systemTB     4.995e-01  8.616e-01  8.979e+00   0.580    0.576
# systemCH     2.395e+00  8.914e-01  8.961e+00   2.687    0.025 *
# seasonYear   4.551e-03  5.704e-03  8.991e+03   0.798    0.425
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#            (Intr) systCK systTB systCH
# systemCK   -0.045
# systemTB   -0.060  0.597
# systemCH   -0.056  0.577  0.689
# seasonYear -0.998  0.000  0.006  0.004

model2 <- lmer(n ~ 1 + system + seasonYear + (1 + seasonYear|systemZone),
               data = modelDF)
summary(model2)






library(broom.mixed)
fit_augmented <- augment(glme_rich_summer)

p1 <-ggplot(fit_augmented, aes(.fitted, .resid))+
  geom_point()

p1 + stat_smooth(method="loess") + 
  geom_hline(yintercept=0, col="red", linetype="dashed") +
  xlab("Fitted values") + 
  ylab("Residuals") + 
  ggtitle("Residual vs Fitted Plot") + 
  theme_bw()

qqnorm(fit_augmented[[".resid"]])

library(aods3)
gof(glme_rich_summer)

sims <- simulate(glme_rich_summer,nsim=1000)
nzeros <- colSums(sims==0)
par(las=1,bty="l")
plot(pt <- prop.table(table(nzeros)),
     ylab="Probability",xlab="Number of zeros")
(obszero <- sum(modelDF$n==0))

#gmod_lme4_agq <- update(gmod_lme4_L,nAGQ=10)



##### brms ####
library(brms)

options(mc.cores = parallel::detectCores())

brm1 <- brm(formula = n ~ system + 
              seasonYear + 
              (1|systemZone) + 
              offset(log(n_hauls)),
            data = modelDF,
            family = "poisson")

(summ_brm1 <- summary(brm1))
plot(brm1)

#####
#define vars for summarizing
vars = c("n", "BottomVegCover")
RichnessSummary <- FullRichness %>%
  subset(select = c(system, seasonYear, season, n, BottomVegCover)) %>%
  group_by(system, seasonYear, season) %>%
  summarise(across(any_of(vars),
                   list(
                     mean = ~mean(., na.rm=TRUE),
                     median = ~median(., na.rm=TRUE),
                     q10 = ~quantile(., probs = 0.1, na.rm=TRUE),
                     q90 = ~quantile(., probs = 0.9, na.rm=TRUE)
                   ), .names = "{col}_{fn}"))

#assign factor levels for plot order (North -> South)
RichnessSummary$system = factor(RichnessSummary$system, levels=c("AP","CK","TB","CH"))

#### Plots ####

library(ggplot2)
library(cowplot)

#Richness Plot
p1 <- ggplot(data=RichnessSummary,
       aes(x=seasonYear, y=n_mean)) +
  geom_ribbon(
    aes(ymin=n_q10,
        ymax=n_q90,
        fill="10th-90th Percentile"),
    linetype=2, alpha=0.1, color="purple") +
  geom_point(
    aes(y = n_median, 
        color="median")
  ) +
  #theme(legend.position="bottom") +
  theme_bw() + 
  ggtitle("Richness over time") +
  xlab("Year") +
  scale_x_continuous(breaks = seq(1998, 2020, 2)) +
  ylab("Number of taxa per haul") +
  geom_line(aes(color = "mean")) +
  facet_grid(season ~ system) +
  theme(strip.text.x = element_text(size = 20)) +
  theme(strip.text.y = element_text(size = 20))

#BottomVegCover plot
p2 <- ggplot(data=RichnessSummary,
       aes(x=seasonYear, y=BottomVegCover_mean)) +
  geom_ribbon(
    aes(ymin=BottomVegCover_q10,
        ymax=BottomVegCover_q90,
        fill="10th-90th Percentile"),
    linetype=2, alpha=0.1, color="purple") +
  geom_point(
    aes(y = BottomVegCover_median, 
        color="median")
  ) +
  #theme(legend.position="bottom") +
  theme_bw() + 
  ggtitle("Seagrass coverage over time") +
  xlab("Year") +
  scale_x_continuous(breaks = seq(1998, 2020, 2)) +
  ylab("Estimated percent seagrass coverage") +
  geom_line(aes(color = "mean")) +
  facet_grid(season ~ system) +
  theme(strip.text.x = element_text(size = 20)) +
  theme(strip.text.y = element_text(size = 20))

plot_grid(p1, p2,
          ncol = 1)
