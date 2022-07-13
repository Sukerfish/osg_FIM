#### Import and Filter Data #####
 
library(tidyverse)
library(lubridate)

load("C:/Users/Gymnothorax/Box/Graduate/RGD/osg_FIM/GearCode20Refresh.Rdata")

#MOI <- c(6:9)

YearFilter <- ZoneFilter %>%
  select(system, StartYear, EndYear) %>%
  unique()

SOI <- c(
"summer",
"winter"
)

HydroList <- TidyHydro %>%
  mutate(month = month(Sampling_Date), year = year(Sampling_Date)) %>%
  left_join(YearFilter) %>%
  mutate(season = if_else(month %in% c(4:5), "spring",
                          if_else(month %in% c(6:9), "summer",
                                  if_else(month %in% c(10:11), "fall",
                                          "winter")))) %>%
  #mutate(seasonYear = if_else(month %in% c(1:3), year - 1, year)) %>%
  # filter(seasonYear >= StartYear - 1) %>%
  # filter(seasonYear != EndYear) %>%
  
  mutate(seasonYear = if_else(month %in% c(12), year + 1, year)) %>%
  filter(season != "winter" | seasonYear != EndYear + 1) %>%
  filter(season != "winter" | seasonYear != StartYear) %>%
  filter(seasonYear >= StartYear) %>%
  filter(seasonYear != EndYear | season != "winter") %>%
  subset(select = -c(StartYear, EndYear)) %>%
  filter(season %in% SOI)

HydroSummary <- HydroList %>%
  #filter(system == "CK") %>%
  group_by(system, seasonYear, season) %>%
  summarise(medTemp    = median(Temperature, na.rm=TRUE),
            meanTemp   = mean(Temperature, na.rm=TRUE),
            lowTemp    = quantile(Temperature, probs = .1, na.rm=TRUE),
            highTemp   = quantile(Temperature, probs = .9, na.rm=TRUE),
            n          = n_distinct(Reference),
  )

HydroSummary$system = factor(HydroSummary$system, levels=c("AP","CK","TB","CH"))

# HydroSummary <- HydroList %>%
#   filter(system == "TB") %>%
#   group_by(system, year) %>%
#   summarise(medTemp = median(Salinity, na.rm=TRUE),
#             meanTemp   = mean(Salinity, na.rm=TRUE),
#             lowTemp    = quantile(Salinity, probs = .1, na.rm=TRUE),
#             highTemp   = quantile(Salinity, probs = .9, na.rm=TRUE),
#             n          = n_distinct(Reference),
#   )

#### Biology ####

# HaulFull <- CleanHRBio %>%
#   inner_join(HydroList) %>%
#   subset(select = -c(Sampling_Date:DissolvedO2)) %>%
#   spread(Scientificname, N2) %>%
#   replace(is.na(.),0)

#### Betadiversity ####

# library(codyn)



#### Plots ####

library(ggplot2)

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





ggplot(data=HydroSummary,
       aes(x=seasonYear, y=meanTemp)) +
  geom_ribbon(
    aes(ymin=lowTemp,
        ymax=highTemp,
        fill="10th-90th Percentile"),
    linetype=2, alpha=0.1, color="purple") +
  geom_point(
    aes(y = medTemp, 
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
