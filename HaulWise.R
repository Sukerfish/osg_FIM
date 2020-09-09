# Pair-wise comparisons for reef populations ----
#
# Author: M. Schram
# Version: 1.1 
# Date: 2020-08-10
#
# This function was designed to work with the species-by-site query in the
# FishSurvey Access database. It was initially written to estimate pair-wise
# comparisons of the original 8 paired reefs.
#
# Spc_Dif <- function(df, opt, par1, par2, comp, SU)
#-----------------------------------------------------
# df   = Data frame where rows are observations and columns are factors and
#           species
#
# opt  = numeric value designating the type of comparison to estimate
#           1: Filters raw data into pairwise subsets (default)
#           2: Estimates pairwise differences in species densities by SU
#           3: Estimates pairwise difference in species richness by SU
#
# par1 = character string corresponding to the column name to be used as
#           the first filtering parameter (e.g. "Distance") 
#
# par2 = character string corresponding to the column name to be used as the 
#           second filtering parameter (e.g. "Proximity") 
#
# comp = character string corresponding to the column name
#           designating values to directly compare (e.g. "Type") 
#
# SU   = character string corresponding to the column name designating the 
#           sampling unit (e.g. "Date")

# PENDING UPDATE TO INCLUDE
# opt  = 4: Estimates pairwise difference in biomass by SU 
#              (NOTE: #4 May require a different function because of new array)
# plot = option to enable/disable generic plot for desired estimate
#           0: No
#           1: Yes (default)
#-----------------------------------------------------
HaulWise <- function(df, 
                     opt = 1, 
                     par1 = "Distance", 
                     par2 = "Proximity", 
                     comp = "Type", 
                     SU = "Date"){
  
  out = list()
  design <- expand(df, df[[par1]], df[[par2]])
  for (i in 1:nrow(design))
  {
    temp_df <- filter(df, df[[par1]] %in% design[i,1] & df[[par2]] %in% design[i,2])
    temp_df <- select(temp_df, which(colSums(temp_df != 0) >0))
    temp <- NULL
    if (opt == 1)
    {
      temp <- temp_df
    }
    else
    {   
      temp_step <- sort(unique(temp_df[[SU]]))
      temp_comp <- sort(unique(temp_df[[comp]]))
      for (j in 1:length(temp_step))
      {
        if(nrow(filter(temp_df, temp_df[[SU]] == temp_step[j] &
                       temp_df[[comp]] == temp_comp[1])) > 0 &
           nrow(filter(temp_df, temp_df[[SU]] == temp_step[j] & 
                       temp_df[[comp]] == temp_comp[2])) > 0)
        {
          if(opt == 2)
          {
            x <- select(filter(temp_df, 
                               temp_df[[SU]] == temp_step[j] &
                                 temp_df[[comp]] == temp_comp[1]),
                        where(is.numeric)) -
              select(filter(temp_df, 
                            temp_df[[SU]] == temp_step[j] &
                              temp_df[[comp]] == temp_comp[2]),
                     where(is.numeric))
            temp <- rbind(temp, x)
          }
          else if (opt == 3)
          {
            temp_art<- filter(temp_df,
                              temp_df[[SU]] == temp_step[j] & 
                                temp_df[[comp]] == temp_comp[1])
            temp_nat<- filter(temp_df, 
                              temp_df[[SU]] == temp_step[j] & 
                                temp_df[[comp]] == temp_comp[2])
            temp_id <- select(temp_art, 
                              c("Year", "Month","Season","Date"))
            temp_art<- rowSums(select(temp_art, 
                                      where(is.numeric)) != 0)
            temp_nat<- rowSums(select(temp_nat, 
                                      where(is.numeric)) != 0)
            x <- cbind(temp_id, "Difference" = temp_art - temp_nat)
            temp <- rbind(temp, x)
          }
          
        }
        else
        {
          temp <- temp
        }
      }
    }    
    if(opt == 2)
    {
      temp <- pivot_longer(temp, everything(), names_to = "Species", values_to = "Difference")    
    }
    else{}
    # [[]] allows accessing list elements w/o $ designation
    # out[[paste(design[i,1], design[i,2], sep = "_")]] <- as.data.frame(temp)
    out[[paste(design[i,1], design[i,2], sep = "_")]] <- temp
  }
  if(opt != 1)
  {
    out[["order"]] <- temp_comp
  }
  return(out)
}
