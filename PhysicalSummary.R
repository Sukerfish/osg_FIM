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

#physical parameters of interest
PaOI <- c(
  "Temperature",
  "Salinity"
  #"pH",
  #"DissolvedO2"
)

#summary statistics of physical parameters of interest
HydroSummary <- HydroList %>%
  group_by(system, seasonYear, season) %>%
  summarise(across(any_of(PaOI),
                   list(
                     mean = ~mean(., na.rm=TRUE),
                     median = ~median(., na.rm=TRUE),
                     q10 = ~quantile(., probs = 0.1, na.rm=TRUE),
                     q90 = ~quantile(., probs = 0.9, na.rm=TRUE)
                   ), .names = "{col}_{fn}"))

#assign factor levels for plot order (North -> South)
HydroSummary$system = factor(HydroSummary$system, levels=c("AP","CK","TB","CH"))

#### Plots ####

library(ggplot2)

#Temperatures Plot
ggplot(data=HydroSummary,
       aes(x=seasonYear, y=Temperature_mean)) +
  geom_ribbon(
    aes(ymin=Temperature_q10,
        ymax=Temperature_q90,
        fill="10th-90th Percentile"),
    linetype=2, alpha=0.1, color="purple") +
  geom_point(
    aes(y = Temperature_median, 
        color="median")
  ) +
  #theme(legend.position="bottom") +
  theme_bw() + 
  ggtitle("Water Temperature Over Time") +
  xlab("Year") +
  scale_x_continuous(breaks = seq(1998, 2020, 2)) +
  ylab("Temperature (°C)") +
  geom_line(aes(color = "mean")) +
  facet_grid(season ~ system) +
  theme(strip.text.x = element_text(size = 20)) +
  theme(strip.text.y = element_text(size = 20))

#Salinity Plot
ggplot(data=HydroSummary,
       aes(x=seasonYear, y=Salinity_mean)) +
  geom_ribbon(
    aes(ymin=Salinity_q10,
        ymax=Salinity_q90,
        fill="10th-90th Percentile"),
    linetype=2, alpha=0.1, color="purple") +
  geom_point(
    aes(y = Salinity_median, 
        color="median")
  ) +
  #theme(legend.position="bottom") +
  theme_bw() + 
  ggtitle("Salinity Over Time") +
  xlab("Year") +
  scale_x_continuous(breaks = seq(1998, 2020, 2)) +
  ylab("Salinity (psu)") +
  geom_line(aes(color = "mean")) +
  facet_grid(season ~ system) +
  theme(strip.text.x = element_text(size = 20)) +
  theme(strip.text.y = element_text(size = 20))

#### old plots ####

# plotHydro <- HydroSummary %>%
#   subset(select = -c(system, n)) %>%
#   gather(key = "stat", value = "temp", -year)
# 
# as.factor(plotHydro$stat)

# ggplot(data=plotHydro,
#        aes(x=year, y=temp, color=stat)) +
#   geom_line() +
#   geom_point()
#   geom_point(aes("meanTemp"))
#   geom_ribbon(
#     aes(ymin="lowTemp",
#         ymax="highTemp"),
#     linetype=2, alpha=0.1, color="purple") +
#   geom_point(
#     aes(y = medTemp)
#   ) +
#   #guides(color = guide_legend) +
#   theme(legend.position="none") +
#   theme_bw() +
#   facet_wrap(as.factor(HydroSummary$system), scales = "free") +
#   ggtitle("Tampa Bay Water Temp") +
#   xlab("Year") +
#   ylab("Temperature (?C)") +
#   geom_line()



