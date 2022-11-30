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
  #mutate(system = factor(system, levels = c("AP", "CK", "TB", "CH"))) %>%
  # mutate(n_hauls = log(n_hauls)) %>%
  # filter(season == "summer") %>%
  # filter(system == "TB") %>%
  filter(season == "winter")

ggplot(modelDF, aes(x=n)) +
  geom_histogram(binwidth=1) +
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


