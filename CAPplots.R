# CAP analyses for each system/season combination

#### data input ######
library(tidyverse)
library(vegan)
library(BiodiversityR)
#library(MASS)

#source ('http://www.davidzeleny.net/anadat-r/doku.php/en:customized_functions:ordicenter?do=export_code&codeblock=0')
ordicenter <- function (ord, groups, display = "sites", w = weights(ord, display), 
                        show.groups, ...) 
{
  weights.default <- function(object, ...) NULL
  pts <- scores(ord, display = display, ...)
  w <- eval(w)
  if (length(w) == 1) 
    w <- rep(1, nrow(pts))
  if (is.null(w)) 
    w <- rep(1, nrow(pts))
  if (!missing(show.groups)) 
  {
    take <- groups %in% show.groups
    pts <- pts[take, , drop = FALSE]
    groups <- groups[take]
    w <- w[take]
  }
  out <- seq(along = groups)
  inds <- names(table(groups))
  for (is in inds) 
  {
    gr <- out[groups == is]
    if (length(gr) > 1)
    {
      X <- pts[gr, ]
      W <- w[gr]
      ave <- apply(X, 2, weighted.mean, w = W)
      vegan:::ordiArgAbsorber(ave[1], ave[2], labels = is, FUN = text, ...)
    }
    if (length(gr) == 1)
    {
      X <- pts[gr, ]
      W <- w[gr]
      vegan:::ordiArgAbsorber(X[1], X[2], labels = is, FUN = text, ...)
    }
  }
  invisible()
}

load('TidyGearCode20.Rdata')

#get the biological data and associated site chars
CleanHauls <- TidyBio %>%
  subset(select = c(Reference, system, season, seasonYear, BottomVegCover, systemZone, Scientificname, N2))
#filter(season == "summer" | season == "winter")
#filter(season == "winter")
#filter(season == "summer")

#get water temp from hydrology dataset
WaterTemp <- HydroList %>%
  subset(select = c(Reference, Temperature))

#collect everything by Reference and spread to long form
SXS_full <- CleanHauls %>%
  mutate(N2 = N2^.25) %>% #fourth-root transform
  group_by(Reference) %>%
  filter(sum(N2)>0) %>% #remove all References with 0 taxa found
  spread(Scientificname,N2) %>%
  ungroup() %>%
  #subset(select = -c(Reference, season, systemZone)) %>%
  replace(is.na(.), 0) %>% #replace all NA values with 0s, i.e. counting as true zero
  group_by(season, system, seasonYear) %>%
  mutate(seasonYear = as.factor(seasonYear)) %>%
  left_join(WaterTemp) %>%
  filter(!is.na(Temperature)) %>% #only excludes 7 sampling events from above pool
  unite("systemSeason",
        system:season,
        sep = "_") %>%
  ungroup()
#summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))


##### setup for CAP loop ####
systemSeason_list <- SXS_full %>%
  select(systemSeason) %>%
  distinct()

CAPsforAll <- list()
for(i in systemSeason_list$systemSeason){
  print(i) #watch progress through list
  
  df <- SXS_full %>% #filter out systemSeason of interest
    filter(systemSeason %in% i)
  
  df_spe <- df %>% #pull out taxa only
    subset(select = -c(systemSeason, seasonYear, Reference, systemZone, BottomVegCover, Temperature)) %>%
    select(which(!colSums(., na.rm=TRUE) %in% 0))
  
  df_env <- data.frame(df %>% #pull out environmental variables
    subset(select = c(systemSeason, seasonYear, BottomVegCover, Temperature)) %>%
    mutate(contYear = as.numeric(as.character(seasonYear))))
  
  bf <- (vegdist(df_spe))^0.5 #Bray-Curtis w/ sqrt to reduce negative eigenvalues
  
  CAP <- CAPdiscrim(bf ~ seasonYear,
                    data = df_env,
                    #add = "lingoes",
                    #parallel = 6,
                    #method="bray", 
                    #permutations = 999
                    )
  
  CAPsforAll[[i]] <- CAP
}


# dbrda <- dbrda(bf ~ Temperature + BottomVegCover,
#                data = df_env)

# example_NMDS=metaMDS(df_spe,k=2,trymax=100)

ordiplot(CAP,
         #display = "species",
         type = "n",
         #xlim = c(-1, 1),
         #ylim = c(-2, 2),
         #axes = TRUE
)
ordiellipse(CAP, df_env$seasonYear,
            #draw = "none",
            col=c(unique(as.numeric(df_env$seasonYear))),
            label = TRUE)



ordisurf(example_NMDS,df_env$Temperature,main="",col="forestgreen")
ordiellipse(CAP, df_env$seasonYear, 
            #draw = "none",
            col=c(unique(as.numeric(df_env$seasonYear))),
            label = TRUE)
ordicenter(example_NMDS,
           groups = df_env$seasonYear,
            #draw = "none",
            col=c(unique(as.numeric(df_env$seasonYear))),
            )

ordiarrows (CAP, 
            groups = df_env$contYear, 
            order.by = df_env$contYear, 
            startmark = 1, label = TRUE, length = .1)
# 
# plot(CAP, display = "sites")

