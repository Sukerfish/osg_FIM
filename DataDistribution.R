# diagnostic histograms for data exploration
#
#
#### data input ######

library(tidyverse)
library(ggplot2)

load('TidyGearCode20.Rdata')
load('SXS_filtered.Rdata')

#### richness ####
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

TotalAbundance <- OnlyFish %>%
  group_by(Reference) %>%
  summarise(n = sum(N2)) %>%
  inner_join(TidyRefsList) %>%
  ungroup() %>%
  bind_rows(NoFish) %>%
  mutate(n = n^.25)

# class(FullRichness$system)
# unique(FullRichness$system)
# as.factor(FullRichness$system)
# 
# , levels = c("AP", "CK", "TB", "CH")
# as.ordered), levels = c("AP", "CK", "TB", "CH")

ggplot(data=FullRichness)+
  geom_histogram(binwidth = 1,
                 aes(x=n)) +
  theme_bw() + 
  ggtitle("Richness per haul by estuary") +
  xlab("Number of taxa in each haul") +
  scale_x_continuous(breaks = seq(0, 100, 5)) +
  ylab("Frequency") +
  facet_grid(rows = vars(FullRichness$season), 
             cols = vars(FullRichness$system)) +
  theme(
    text=element_text(size=20))
        
#total abundance
ggplot(data=TotalAbundance)+
  geom_histogram(binwidth = 1,
                 aes(x=n)) +
  theme_bw() + 
  ggtitle("Total abundance per haul by estuary") +
  xlab("Number of individuals in each haul (fourth-root transformed)") +
  scale_x_continuous(breaks = seq(0, 100, 3)) +
  ylab("Frequency") +
  facet_grid(rows = vars(TotalAbundance$season), 
             cols = vars(TotalAbundance$system)) +
  theme(
    text=element_text(size=20))

# p + theme(text=element_text(size=20), #change font size of all text
#           axis.text=element_text(size=20), #change font size of axis text
#           axis.title=element_text(size=20), #change font size of axis titles
#           plot.title=element_text(size=20), #change font size of plot title
#           legend.text=element_text(size=20), #change font size of legend text
#           legend.title=element_text(size=20)) #change font size of legend title 
# 

stuff <- FullRichness %>%
  filter(system == "CK") %>%
  filter(season == "summer")

library(EnvStats)

gof.list <- gofTest(stuff$n)
plot(gof.list)

