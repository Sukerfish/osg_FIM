# Convert NODC Code to Human Readable ----
#
# Author: G. Miller
# Version: 1.0 
# Date: 2020-09-25
#
# This functions works with FWRI FIM data to combine taxa that are difficult
# to ID but can be merged with others one taxonomic designation up. It can also
# convert to human readable names.
#
# osg_CleanBio <- function(df, HRSpecies)
#-----------------------------------------------------
# df    = Data frame of entire hsdb_tbl_corp_biology_number table
#
# HRSpecies = optional list of species names and NODCCODEs
#             REQUIRES column named NODCCODE
#-----------------------------------------------------
osg_ComBio <- function(df,
                       HRSpecies){
  
  HRSpecies$NODCCODE <- as.character(HRSpecies$NODCCODE)
  
  # NODC Codes to combine based on
  # Difficult-to-ID taxa (DTI)
  # Anchoa spp. is        8747020200
  # Eucinostomous spp. is 8835390100

  AncSpp <- c(8747020209, 8747020203, 8747020201, 8747020204, 8747020205, 8747020202, 8747020200)
  EucSpp <- c(8835390101, 8835390102, 8835390111, 8835390108, 8835390109, 8835390104, 8835390105, 8835390100)

  # replace all DTI txa with string then replace that string with the appropriate NODC Code
  # then group by reference and NODC Code to combine duplicates
  # i.e., combines multiple reports of same NODC Code per reference (haul)
  # covers modified codes above and multiple counts of taxa (as with LaRh)
  temp_df <- df %>%
    mutate(NODCCODE=replace(NODCCODE, NODCCODE %in% AncSpp, "AncSpp")) %>%
    mutate(NODCCODE=replace(NODCCODE, NODCCODE %in% EucSpp, "EucSpp")) %>%
    mutate(NODCCODE=replace(NODCCODE, NODCCODE == "AncSpp", 8747020200)) %>%
    mutate(NODCCODE=replace(NODCCODE, NODCCODE == "EucSpp", 8835390100)) %>%
    group_by(Reference, NODCCODE) %>%
    mutate(N2 = sum(N)) %>%
    select(-N) %>%
    distinct() %>%
    ungroup()

  if(missing(HRSpecies)){
    return(temp_df)
  }
  # final lookup of NODC Codes against human readable names
  out <- inner_join(temp_df, HRSpecies, by = "NODCCODE") %>%
    select(Reference, Scientificname, month, year, Stratum, Zone,
           Grid, BottomVegCover,BycatchQuantity, Bank, 
           ShoreDistance, N2)
  return(out)
  }