#### Import and Filter Data #####
 
library(tidyverse)
library(lubridate)

#load("C:/Users/Gymnothorax/Box/Graduate/RGD/osg_FIM/GearCode20.Rdata")
load("C:/Users/Gymnothorax/Box/Graduate/RGD/osg_FIM/GearCode2023160300.Rdata")

MOI <- c(6:9)

YearFilter <- ZoneFilter %>%
  select(system, StartYear, EndYear) %>%
  unique()

SOI <- c(
 # "summer"
  "winter"
 )

HydroSummary <- TidyHydro %>%
  mutate(month = month(Sampling_Date), year = year(Sampling_Date)) %>%
  left_join(YearFilter) %>%
  mutate(season = if_else(month %in% c(4:5), "spring",
                          if_else(month %in% c(6:9), "summer",
                                  if_else(month %in% c(10:11), "fall",
                                          "winter")))) %>%
  mutate(seasonYear = if_else(month %in% c(1:3), year - 1, year)) %>%
  #filter(seasonYear >= StartYear - 1) %>%
  #filter(seasonYear != EndYear) %>%
  subset(select = -c(StartYear, EndYear)) %>%
  #filter(season %in% SOI) %>%
  filter(system == "TB") %>%
  group_by(system, year) %>%
  summarise(medTemp = median(Temperature, na.rm=TRUE),
            meanTemp   = mean(Temperature, na.rm=TRUE),
            lowTemp    = quantile(Temperature, probs = .1, na.rm=TRUE),
            highTemp   = quantile(Temperature, probs = .9, na.rm=TRUE),
            n          = n_distinct(Reference),
  )

#### Plots ####

library(ggplot2)

plotHydro <- HydroSummary %>%
  subset(select = -c(system, n)) %>%
  gather(key = "stat", value = "temp", -year)

as.factor(plotHydro$stat)

ggplot(data=plotHydro,
       aes(x=year, y=temp, color=stat)) +
  geom_line() +
  geom_point()
  geom_point(aes("meanTemp"))
  geom_ribbon(
    aes(ymin="lowTemp",
        ymax="highTemp"),
    linetype=2, alpha=0.1, color="purple") +
  geom_point(
    aes(y = medTemp)
  ) +
  #guides(color = guide_legend) +
  theme(legend.position="none") +
  theme_bw() + 
  facet_wrap(as.factor(HydroSummary$system), scales = "free") +
  ggtitle("Tampa Bay Water Temp") +
  xlab("Year") +
  ylab("Temperature (°C)") +
  geom_line()





ggplot(data=HydroSummary,
       aes(x=year, y=meanTemp)) +
  geom_ribbon(
    aes(ymin=lowTemp,
        ymax=highTemp),
    linetype=2, alpha=0.1, color="purple") +
  geom_point(
    aes(y = medTemp)
  ) +
  #guides(color = guide_legend) +
  theme(legend.position="none") +
  theme_bw() + 
  #facet_wrap(as.factor(HydroSummary$system), scales = "free") +
  ggtitle("Tampa Bay Water Temp") +
  xlab("Year") +
  ylab("Temperature (°C)") +
  geom_line()
