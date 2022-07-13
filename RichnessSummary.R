#### Import and Filter Data #####
 
library(tidyverse)
library(lubridate)

load('GearCode20Refresh.Rdata')

#establish years of interest
YearFilter <- ZoneFilter %>%
  select(system, StartYear, EndYear) %>%
  unique()

#establish seasons of interest
SOI <- c(
  "summer",
  "winter"
)

#expand sampling date and attach years of interest
HydroList <- TidyHydro %>%
  mutate(month = month(Sampling_Date), year = year(Sampling_Date)) %>%
  left_join(YearFilter) %>%
  #label specific months as seasons
  mutate(season = if_else(month %in% c(4:5), "spring",
                          if_else(month %in% c(6:9), "summer",
                                  if_else(month %in% c(10:11), "fall",
                                          "winter")))) %>%
  #tag everything with year associated with sampling season, with December
  #as the preceding year's data
  mutate(seasonYear = if_else(month %in% c(12), year + 1, year)) %>%
  #remove the incomplete final seasonYear of winter data
  filter(season != "winter" | seasonYear != EndYear + 1) %>%
  #filter by years of interest (all begin before, so first year of interest is
  #not truncated)
  filter(seasonYear >= StartYear) %>%
  #remove year of interest logic
  subset(select = -c(StartYear, EndYear)) %>%
  #filter by seasons of interest
  filter(season %in% SOI) 

TidyRefsList <- HydroList %>%
  subset(select = c(Reference, system, season, seasonYear))

TidyBio <- TidyRefsList %>%
  left_join(CleanHRBio, by = "Reference")

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
  ungroup %>%
  bind_rows(NoFish)
  
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


