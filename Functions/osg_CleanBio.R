# Pair-wise comparisons for reef populations ----
#
# Author: G. Miller
# Version: 1.0 
# Date: 2020-09-25
#
# This functions works with FWRI FIM data to remove a few problematic taxa
# and to select samples of interest for further analyses
#
# osg_CleanBio <- function(df, opt, par1, par2, comp, SU)
#-----------------------------------------------------
# df    = Data frame of entire hsdb_tbl_corp_biology_number table
#
# RefOI = list of sampling units of interest as Reference numbers
#-----------------------------------------------------
osg_CleanBio <- function(df, 
                     RefOI){
  
  out = data.frame()
  out <- df %>%
    select(Reference, Species_record_id, Splittype, Splitlevel, Cells,
           NODCCODE, Number, FHC) %>%
    filter(FHC != "D" &
           #remove unidentified species
           !NODCCODE %in% c("9999000000") &
           #remove freshwater prawns
           !NODCCODE %in% c("6179110200", "6179110201")) %>%
    collect() %>%
    filter(#remove all turtles, by NODCCODEs that start with 9002
    !str_detect(NODCCODE, "^9002"),
    #remove all grass shrimp, by NODCCODEs that start with 617911
    !str_detect(NODCCODE, "^617911")) %>%
    #subset to references of interest
    inner_join(RefOI, by = "Reference") %>%
    #account for splitter
    mutate(Count = case_when(!is.na(as.numeric(Splittype)) ~ Number*(as.numeric(Splittype)^as.numeric(Splitlevel)),
                           TRUE ~ as.numeric(Number))) %>%
    select(-Splittype, -Splitlevel, -Species_record_id, -Cells, -FHC) %>%
    group_by(Reference, NODCCODE) %>%
    mutate(N = sum(Count)) %>%
    select(-Count, -Number) %>%
    distinct() %>%
    ungroup() %>%
    #select only a few variables and expand sampling date
    select(Reference, NODCCODE, Sampling_Date, Stratum, Zone,
           Grid, BottomVegCover,BycatchQuantity, Bank, ShoreDistance, N) %>%
    mutate(month = month(Sampling_Date), year = year(Sampling_Date)) %>%
    select(-Sampling_Date)
  
  return(out)
  }
