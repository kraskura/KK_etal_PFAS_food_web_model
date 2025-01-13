
library(tidyverse)
library(readxl)
library(here)
here::i_am(path = "./Code/R data tables/WG_fallWinter.R")

# seasons 
# JBA/AFBC Spring sampling (water data): March 30  - April 28 2021: 
# JBA/AFBC Summer sampling (water data): June 14 - July 16 2021
# JBA/AFBC Fall/Winter sampling (water data): November 19 2020 - January 14 2021

# WG/NRP Fall sampling (water data): October 9 - November 10 2020; second to last day 2020-11-05
# WG/NRP Summer sampling (water data): June 15 - August 19 (some sparse sampling) 2021
# WG/NRP Spring sampling (water data): February 9 - March 23 2021

# sample dates
# WG/NRP   as.Date("2020-11-10") (fall/winter) # <<< this file, median values
            # as.Date("2021-04-28") (spring)
            # as.Date("2021-08-19") (summer) 
 
# JBA/AFBC as.Date("2021-04-28") (spring)
         # as.Date("2021-07-16") (summer)

# pull in relevant data:
d.fish<-read_xlsx(here("Data", "FINAL", "WG_Biota_Data.xlsx"), sheet = "Data") #fish
d.water<-read_xlsx(here("Data", "FINAL", "WG_Water_Data.xlsx"), sheet = "Data") #water
d.sed<-read_xlsx(here("Data", "FINAL", "WG_Sediment_Data.xlsx"), sheet = "Data") #sediment

# filter out data needed and format. 
d.fish<-d.fish %>%
  mutate_at(vars(Species, Tissue), factor) %>%
  mutate(Seasonality = if_else(SampleDate > as.Date("2021-06-01"), "summer", 
                              if_else(SampleDate < as.Date("2020-11-20"), "fall/winter",
                                      "spring"))) %>% 
  filter(c(Tissue == "Whole Body" | Tissue == "Muscle") & Seasonality == "fall/winter") %>% 
  select(Species, Tissue, Weight, 
         PFBA, PFPeA, PFHxA, PFHpA, PFOA, PFOS,
         PFNA, PFHpA, PFHxA, PFHxS, PFPeS) %>% 
  mutate(PFUA = 0, PFDA = 0) %>% 
  mutate(Sp = if_else(Species == "Banded Killifish", "Kil",
                if_else(Species == "Creek Chubsucker", "Chu", 
                   if_else(Species == "Dace sp.", "Dac",
                      if_else(Species == "Darter sp.", "Dar",
                        if_else(Species == "Eastern Mudminnow", "Min",
                          if_else(Species == "Margined Madtom", "Mad",
                            if_else(Species == "Pumpkinseed", "Pum",
                              if_else(Species == "Swallowtail Shiner", "Swa",
                                if_else(Species == "Fallfish","Fal",
                                  if_else(Species == "Largemouth Bass", "Bas",
                                    if_else(Species == "Bluegill", "Bgl",
                                      if_else(Species == "Prey", "Pry", NA)))))))))))))

# Fish masses
d.fish.w <- d.fish %>% 
  group_by(Sp, Species) %>% 
  summarise(n = n(), 
            median_mass = median(Weight, na.rm = TRUE), 
            min_mass = min(Weight, na.rm = TRUE), 
            max_mass = max(Weight, na.rm = TRUE))

# fillet-to-whole body conversion 
d.fish[d.fish$Tissue == "Muscle", c(4:14)]<-d.fish[d.fish$Tissue == "Muscle", c(4:14)]*2.5
d.fish<-d.fish %>% 
  dplyr::group_by(Sp) %>% 
  summarize(PFHxS.m = median(PFHxS), PFHxS.min = min(PFHxS), PFHxS.max = max(PFHxS),
            PFOS.m = median(PFOS), PFOS.min = min(PFOS), PFOS.max = max(PFOS),
            PFOA.m = median(PFOA), PFOA.min = min(PFOA), PFOA.max = max(PFOA),
            PFNA.m = median(PFNA), PFNA.min = min(PFNA), PFNA.max = max(PFNA),
            PFDA.m = median(PFDA), PFDA.min = min(PFDA), PFDA.max = max(PFDA),
            PFUA.m = median(PFUA), PFUA.min = min(PFUA), PFUA.max = max(PFUA), 
            n = n())

d.water<-d.water %>%
  mutate_at(vars(SampleID, Location), factor) %>%
  mutate(Seasonality = if_else(SampleDate > as.Date("2021-06-01"), "summer", 
                              if_else(SampleDate < as.Date("2020-11-20"), "fall/winter",
                                      "spring"))) %>% 
  filter(Seasonality == "fall/winter") %>% 
  select(SampleDate, Location,
         WaterTemp,
         PFBA, PFPeA, PFHxA, PFHpA, PFOA, PFOS,
         PFNA, PFHpA, PFHxA, PFHxS, PFPeS) %>%
  mutate(PFUA = 0, PFDA = 0) %>% 
  summarize(PFHxS.m = median(PFHxS), PFHxS.min = min(PFHxS), PFHxS.max = max(PFHxS),
            PFOS.m = median(PFOS), PFOS.min = min(PFOS), PFOS.max = max(PFOS),
            PFOA.m = median(PFOA), PFOA.min = min(PFOA), PFOA.max = max(PFOA),
            PFNA.m = median(PFNA), PFNA.min = min(PFNA), PFNA.max = max(PFNA),
            PFDA.m = median(PFDA), PFDA.min = min(PFDA), PFDA.max = max(PFDA),
            PFUA.m = median(PFUA), PFUA.min = min(PFUA), PFUA.max = max(PFUA),
            temp.m = median(WaterTemp), temp.min = min(WaterTemp), temp.max = max(WaterTemp), 
            n = n())

d.sed<-d.sed%>%
  mutate_at(vars(SampleID, Location), factor) %>% 
  mutate(Seasonality = if_else(SampleDate > as.Date("2021-06-01"), "summer", 
                              if_else(SampleDate < as.Date("2020-11-20"), "fall/winter",
                                      "spring"))) %>% 
  filter(Seasonality == "fall/winter") %>% 
  select(SampleDate, Location,
         PFBA, PFPeA, PFHxA, PFHpA, PFOA, PFOS,
         PFNA, PFHpA, PFHxA, PFHxS, PFPeS) %>% 
  mutate(PFUA = 0, PFDA = 0) %>% 
  summarize(PFHxS.m = median(PFHxS), PFHxS.min = min(PFHxS), PFHxS.max = max(PFHxS),
            PFOS.m = median(PFOS), PFOS.min = min(PFOS), PFOS.max = max(PFOS),
            PFOA.m = median(PFOA), PFOA.min = min(PFOA), PFOA.max = max(PFOA),
            PFNA.m = median(PFNA), PFNA.min = min(PFNA), PFNA.max = max(PFNA),
            PFDA.m = median(PFDA), PFDA.min = min(PFDA), PFDA.max = max(PFDA),
            PFUA.m = median(PFUA), PFUA.min = min(PFUA), PFUA.max = max(PFUA), 
            n = n())



# **********************************
# **********************************
# order for PFAA is c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")
inputFiles_list<-create_data_tables(
  species = c("Phy", "Pry", "Bgl", "Bas"), 
  group_species = c("plant", "fish", "fish", "fish"), 
  WB_kg = unlist(c(NA, 
            5/1000, # prey = 5 g (convert in kg)
            d.fish.w[d.fish.w$Sp == "Bgl", "median_mass"]/1000,
            d.fish.w[d.fish.w$Sp == "Bas", "median_mass"]/1000)), 
  m_O = c(1, 1, 1, 1),
  GRF = c(0.8, 0.00150, 0.00150, 0.00150),
  P_B = c(0.5, 0.15, 0.15, 0.15), 
  diet = list(
    Pry	= unlist(c(d.fish[d.fish$Sp == "Pry", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Pry", "PFOS.m"],
            d.fish[d.fish$Sp == "Pry", "PFOA.m"],
            d.fish[d.fish$Sp == "Pry", "PFNA.m"],
            d.fish[d.fish$Sp == "Pry", "PFDA.m"],
            d.fish[d.fish$Sp == "Pry", "PFUA.m"])),
    Bgl	= unlist(c(d.fish[d.fish$Sp == "Bgl", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Bgl", "PFOS.m"],
            d.fish[d.fish$Sp == "Bgl", "PFOA.m"],
            d.fish[d.fish$Sp == "Bgl", "PFNA.m"],
            d.fish[d.fish$Sp == "Bgl", "PFDA.m"],
            d.fish[d.fish$Sp == "Bgl", "PFUA.m"])), # same
    Bas	= unlist(c(d.fish[d.fish$Sp == "Bas", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Bas", "PFOS.m"],
            d.fish[d.fish$Sp == "Bas", "PFOA.m"],
            d.fish[d.fish$Sp == "Bas", "PFNA.m"],
            d.fish[d.fish$Sp == "Bas", "PFDA.m"],
            d.fish[d.fish$Sp == "Bas", "PFUA.m"]))),
  min_diet = list(
    Pry	= unlist(c(d.fish[d.fish$Sp == "Pry", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Pry", "PFOS.min"],
            d.fish[d.fish$Sp == "Pry", "PFOA.min"],
            d.fish[d.fish$Sp == "Pry", "PFNA.min"],
            d.fish[d.fish$Sp == "Pry", "PFDA.min"],
            d.fish[d.fish$Sp == "Pry", "PFUA.min"])),
    Bgl	= unlist(c(d.fish[d.fish$Sp == "Bgl", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Bgl", "PFOS.min"],
            d.fish[d.fish$Sp == "Bgl", "PFOA.min"],
            d.fish[d.fish$Sp == "Bgl", "PFNA.min"],
            d.fish[d.fish$Sp == "Bgl", "PFDA.min"],
            d.fish[d.fish$Sp == "Bgl", "PFUA.min"])), # same
    Bas	= unlist(c(d.fish[d.fish$Sp == "Bas", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Bas", "PFOS.min"],
            d.fish[d.fish$Sp == "Bas", "PFOA.min"],
            d.fish[d.fish$Sp == "Bas", "PFNA.min"],
            d.fish[d.fish$Sp == "Bas", "PFDA.min"],
            d.fish[d.fish$Sp == "Bas", "PFUA.min"]))),
  max_diet = list(
    Pry	= unlist(c(d.fish[d.fish$Sp == "Pry", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Pry", "PFOS.max"],
            d.fish[d.fish$Sp == "Pry", "PFOA.max"],
            d.fish[d.fish$Sp == "Pry", "PFNA.max"],
            d.fish[d.fish$Sp == "Pry", "PFDA.max"],
            d.fish[d.fish$Sp == "Pry", "PFUA.max"])),
    Bgl	= unlist(c(d.fish[d.fish$Sp == "Bgl", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Bgl", "PFOS.max"],
            d.fish[d.fish$Sp == "Bgl", "PFOA.max"],
            d.fish[d.fish$Sp == "Bgl", "PFNA.max"],
            d.fish[d.fish$Sp == "Bgl", "PFDA.max"],
            d.fish[d.fish$Sp == "Bgl", "PFUA.max"])), # same
    Bas	= unlist(c(d.fish[d.fish$Sp == "Bas", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Bas", "PFOS.max"],
            d.fish[d.fish$Sp == "Bas", "PFOA.max"],
            d.fish[d.fish$Sp == "Bas", "PFNA.max"],
            d.fish[d.fish$Sp == "Bas", "PFDA.max"],
            d.fish[d.fish$Sp == "Bas", "PFUA.max"]))),
  foodWeb = list(
Phy = c(1,0,0,0,0),
Pry = c(0.7,0.3,0,0,0),
Bgl = c(0.1,0.3,0.6,0,0),
Bas = c(0.1,0.1,0.7,0.1,0)
      ),
  C_WTO_ng_mL = unlist(c(d.water[, "PFHxS.m"]/1000, # original data in ng/L convert to ng/mL
                  d.water[, "PFOS.m"]/1000,
                  d.water[, "PFOA.m"]/1000,
                  d.water[, "PFNA.m"]/1000,
                  d.water[, "PFDA.m"]/1000,
                  d.water[, "PFUA.m"]/1000)), 
  C_WTO_max_ng_mL = unlist(c(d.water[, "PFHxS.max"]/1000, 
                  d.water[, "PFOS.max"]/1000,
                  d.water[, "PFOA.max"]/1000,
                  d.water[, "PFNA.max"]/1000,
                  d.water[, "PFDA.max"]/1000,
                  d.water[, "PFUA.max"]/1000)), 
  C_WTO_min_ng_mL = unlist(c(d.water[, "PFHxS.min"]/1000, 
                  d.water[, "PFOS.min"]/1000,
                  d.water[, "PFOA.min"]/1000,
                  d.water[, "PFNA.min"]/1000,
                  d.water[, "PFDA.min"]/1000,
                  d.water[, "PFUA.min"]/1000)), 
  C_s_ng_g = unlist(c(d.sed[, "PFHxS.m"], # in ng/g
                  d.sed[, "PFOS.m"],
                  d.sed[, "PFOA.m"],
                  d.sed[, "PFNA.m"],
                  d.sed[, "PFDA.m"],
                  d.sed[, "PFUA.m"])), # order of PFAA
  C_s_max_ng_g = unlist(c(d.sed[, "PFHxS.max"], 
                  d.sed[, "PFOS.max"],
                  d.sed[, "PFOA.max"],
                  d.sed[, "PFNA.max"],
                  d.sed[, "PFDA.max"],
                  d.sed[, "PFUA.max"])), # order of PFAA
  C_s_min_ng_g = unlist(c(d.sed[, "PFHxS.min"], 
                  d.sed[, "PFOS.min"],
                  d.sed[, "PFOA.min"],
                  d.sed[, "PFNA.min"],
                  d.sed[, "PFDA.min"],
                  d.sed[, "PFUA.min"])),# order of PFAA
  C_OX = 8, # assumed
  T = unlist(c(d.water[, "temp.m"])),
  # c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")
  log_Koc = c("Brown etal", 2.3317871	, 2.891004, 2.3032831,	3.0329491,	6,	5) # c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")

  # log_Dmw = c("calibr-cosmotherm", 3.37, 4.61, 3.47, 4.10, 4.69, 5.34)
           			
) 


