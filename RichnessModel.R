# #### Import and Filter Data #####

library(tidyverse)

load('TidyGearCode20.Rdata')

OnlyFish <- TidyBio %>%
  filter(Scientificname != "No fish")

NoFish <- TidyBio %>%
  filter(Scientificname == "No fish") %>%
  mutate(n = 0) %>%
  subset(select = c(Reference, n, system, season, seasonYear))

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

modelDF <- RefsList %>%
  subset(select = c(Reference, Zone)) %>%
  inner_join(FullRichness) %>%
  group_by(seasonYear, system, season, Zone) %>%
  add_count(name = "n_hauls") %>%
  ungroup() %>%
  mutate(systemZone = str_c(system, "_", Zone)) %>%
  filter(season == "summer")

unique(modelDF$systemZone)
unique(modelDF$Zone)
str_to_upper(modelDF$systemZone, locale = "en")

library(nlme)
library(MASS)

glmmPQL1 <- glmmPQL(fixed = n ~ system + seasonYear + offset(log(n_hauls)),
                    random = ~ 1|systemZone,
                    family = "poisson",
                    correlation=corARMA(form=~1|systemZone/seasonYear,p=1),
                    data = modelDF)
summary(glmmPQL1)
plot(glmmPQL1)

#kind of works, cant converge
gmod_lme4_L <- glmer(n~factor(system)+
                       offset(log(n_hauls))+
                       seasonYear+
                       (1|systemZone),
                     family=poisson,
                     data=modelDF,
                     control=glmerControl(optimizer="bobyqa",
                                          check.conv.grad=.makeCC("warning",0.05)))
summary(gmod_lme4_L)
gmod_lme4_agq <- update(gmod_lme4_L,nAGQ=10)

#####
options(mc.cores = parallel::detectCores())
library(brms)
brm1 <- brm(formula = n ~ system + 
              seasonYear + 
              (1|Zone) + 
              offset(log(n_hauls)),
            data = modelDF,
            family = "poisson")

(summ_brm1 <- summary(brm1))
plot(brm1)

#####
library(lme4)

m2 <- glm(modelDF$n ~ modelDF$seasonYear + offset(log(offsetDF)), family = poisson, data = df)


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


