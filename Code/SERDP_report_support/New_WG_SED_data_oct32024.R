# library(ggplot2)
library(readxl)
library(tidyverse)
library(here)

# existing data:
d_wg_sed<-read.csv(here("./Data/SERDP_report_support/WG_Sediment_Data.csv"))

# new data: 
d0<-read_excel(here("./Data/SERDP_report_support/SGS PFAS includes some SpringWG sediment data.xlsx"), sheet = 1)
unique(d0$SAMPLE_NO)
removals<-unique(d_wg_sed$SampleID)
# the ones currently 

d0 %>% 
  filter(UNIT != "% Recovery" &
           SAMPLE_NO != "Lab Blank" &
           UNIT != "%") %>% 
  select(SAMPLE_NO, COMPOUND, CONC_FOUND, DETECTION_LIMIT) %>% 
  dplyr::rename("SampleID" = "SAMPLE_NO") %>% 
  filter(!SampleID %in% removals) %>% 
  filter(grepl("WG", SampleID)) %>% 
  mutate("Location" = gsub(replacement = "", "-", x = substr(SampleID, start = 11, stop = 15), ),
         "SampleDate" = as.Date(
           substr(SampleID, nchar(SampleID) +1 - 8, nchar(SampleID)),
           format = "%m%d%Y"), 
         "SedimentType" = NA, 
         "GrainSize" = NA,
         "TOC" = NA,  
         "Gravel" = NA, 
         "Silt_Clay" = NA, 
         "Solids" = NA,
         "SumPFAS" = NA) %>% 
  pivot_wider(values_from = CONC_FOUND,
              names_from = COMPOUND)
  
         


# all new samples have detection limits < identified compound concentrations
#  [1] "SampleID"         "Location"         "SampleDate"       "SedimentType"     "GrainSize"        "TOC"             
#  [7] "Gravel"           "Silt_Clay"        "Solids"           "PFBA"             "PFPeA"            "PFBS"            
# [13] "PFHxA"            "PFPeS"            "PFHpA"            "PFHxS"            "PFOA"             "X6.2.FTS"        
# [19] "PFHpS"            "PFNA"             "PFOS"             "PFDA"             "PFUnA"            "PFDoA"           
# [25] "PFTrDA"           "PFTeDA"           "PFNS"             "PFDoS"            "X4.2.FTS"         "X8.2.FTS"        
# [31] "FOSA"             "N.MeFOSA"         "N.MeFOSAA"        "N.EtFOSA"         "N.EtFOSAA"        "N.MeFOSE"        
# [37] "N.EtFOSE"         "HFPO.DA"          "ADONA"            "PFMBA"            "NFDHA"            "X9Cl.PF3ONS"     
# [43] "X11Cl.PF3OUdS"    "PFEESA"           "X3.3.FTCA"        "X5.3.FTCA"        "X7.3.FTCA"        "PFDS"            
# [49] "PFMPA"            "SumPFAS"