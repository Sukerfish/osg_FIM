# diagnostic histograms for data exploration
#
#
#### data input ######

library(tidyverse)
library(ggplot2)

load('TidyGearCode20.Rdata')

#### richness ####
CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, systemZone, Scientificname, N2))
  #filter(season == "summer" | season == "winter")
 # filter(season == "winter")
#filter(season == "summer")

###### main ######
SiteXSpeciesFull <- CleanHauls %>%
  mutate(N2 = N2^.25) %>% #fourth-root transform
  group_by(Reference) %>%
  spread(Scientificname,N2) %>%
  ungroup() %>%
  subset(select = -c(Reference, season, systemZone)) %>%
  replace(is.na(.), 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  group_by(system, seasonYear) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))


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

