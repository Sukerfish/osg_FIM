# #### Import and Filter Data #####

library(tidyverse)

load('TidyGearCode20.Rdata')

OnlyFish <- TidyBio %>%
  filter(Scientificname != "No fish")

NoFish <- TidyBio %>%
  filter(Scientificname == "No fish") %>%
  mutate(n = 0) %>%
  subset(select = c(Reference, n, system, season, seasonYear))

FullAb <- OnlyFish %>%
  group_by(Reference) %>%
  summarise(n = sum(N2)) %>%
  inner_join(TidyRefsList) %>%
  ungroup() %>%
  bind_rows(NoFish) %>%
  mutate(n = n^.25)
  
vars = c("n")

AbSummary <- FullAb %>%
  group_by(system, seasonYear, season) %>%
  summarise(across(any_of(vars),
                   list(
                     mean = ~mean(., na.rm=TRUE),
                     median = ~median(., na.rm=TRUE),
                     q10 = ~quantile(., probs = 0.1, na.rm=TRUE),
                     q90 = ~quantile(., probs = 0.9, na.rm=TRUE)
                   ), .names = "{col}_{fn}"))

#assign factor levels for plot order (North -> South)
AbSummary$system = factor(AbSummary$system, levels=c("AP","CK","TB","CH"))

#### Plots ####

library(ggplot2)

#TotalAb Plot
AbSummaryPlot <- ggplot(data=AbSummary,
       aes(x=seasonYear, y=n_mean)) +
  geom_ribbon(
    aes(ymin=n_q10,
        ymax=n_q90,
        fill="10th-90th percentile"),
    linetype=2, alpha=0.1, color="purple") +
  geom_point(
    aes(y = n_median, 
        color="median")
  ) +
  geom_line(aes(color = "mean")) +
  theme_bw() + 
  labs(
    title = "Total abundance over time",
    x = "Year",
    y = "Number of individuals per haul \n (4th root transformed)",
  ) +
  scale_x_continuous(breaks = seq(1998, 2020, 5)) +
  facet_grid(season ~ system, scales = "free_x") +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    #strip.text.x = element_text(size = 20),
    #strip.text.y = element_text(size = 20),
    plot.title = element_text(face="bold"),
    text = element_text(family="serif"))


