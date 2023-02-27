library(tidyverse)

load('TidyGearCode20.Rdata')

taxaList <- TidyBio %>%
  select(Scientificname) %>%
  distinct()

write.csv(taxaList, './Outputs/taxaList.csv')
