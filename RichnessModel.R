# #### Import and Filter Data #####

library(tidyverse)

load('TidyGearCode20.Rdata')

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
 
# GLMM:
# Response: richness or total abundance
# Fixed effects: estuary (four levels) and year
# Random effect: within-estuary zone
# Offset term for variation in number of seine hauls
# Autoregressive term to account for serial correlation

modelDF <- FullRichness %>%
  group_by(seasonYear, system, season) %>%
  add_count(name = "n_hauls") %>%
  ungroup() %>%
  mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH"))) %>%
  # mutate(n_hauls = log(n_hauls)) %>%
  # filter(season == "summer") %>%
  # filter(system == "TB") %>%
  filter(season == "winter")

ggplot(modelDF, aes(x=n)) +
  geom_histogram(binwidth=1) +
  theme(axis.text=element_text(size = 12)) +
  theme(axis.title=element_text(size = 16)) +
  theme(strip.text = element_text(size = 16)) +
  ggtitle("Richness per haul") +
  theme(title=element_text(size = 20)) +
  facet_grid(season~system)

ggplot(modelDF, aes(y=n,
                    x = factor(seasonYear))) +
  geom_boxplot() +
  scale_x_discrete(breaks=c(2000,2005,2010,2015,2020))+
  theme(axis.text=element_text(size = 12)) +
  theme(axis.title=element_text(size = 16)) +
  theme(strip.text = element_text(size = 16)) +
  ggtitle("Richness per haul over time") +
  theme(title=element_text(size = 20)) +
  facet_grid(season~system)

modelDF$systemZone <- as.factor(modelDF$systemZone)

# test <-modelDF %>%
#   subset(select = c(seasonYear))
# df <- fct_reorder(factor(modelDF$seasonYear),test$seasonYear)

###### pql ####
library(nlme)
library(lme4)
library(MASS) # needs MASS (version 7.3-58)

#DOES NOT WORK
glmmPQL <- glmmPQL(n ~ system + seasonYear + offset(log(n_hauls)),
                    random = ~ 1|systemZone,
                    family = poisson,
                    correlation = corARMA(form = ~1|systemZone, p = 1, q = 1),
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
rmod_tmb <- glmmTMB(n~system+
                      seasonYear+
                      offset(log(n_hauls))+
                      ar1(factor(seasonYear) + 1|systemZone)+
                      (1|systemZone),
                    zi=~0,
                    family=poisson,
                    data=modelDF)
summary(rmod_tmb)
VarCorr(rmod_tmb)

##### lme4 ####
#converges if seasonYear is random
library(lme4)
library(lmerTest)

gmod_lme4_L <- glmer(n~system+
                       offset(log(n_hauls))+
                       (1|systemZone/seasonYear),
                     family="poisson",
                     data=modelDF,
                     control=glmerControl(optimizer="bobyqa",
                                          check.conv.grad=.makeCC("warning",0.05)))
summary(gmod_lme4_L)
plot(gmod_lme4_L)
plot(gmod_lme4_L,systemZone~resid(.))

gmod_lme4_L <- glmer(n~system+seasonYear+
                       offset(log(n_hauls))+
                       (1|systemZone),
                     family="poisson",
                     data=modelDF,
                     control=glmerControl(optimizer="bobyqa",
                                          check.conv.grad=.makeCC("warning",0.05)))
summary(gmod_lme4_L)


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
fit_augmented <- augment(gmod_lme4_L)

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
gof(gmod_lme4_L)

sims <- simulate(gmod_lme4_L,nsim=1000)
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

vars = c("n")

RichnessSummary <- FullRichness %>%
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

#Richness Plot
ggplot(data=RichnessSummary,
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


