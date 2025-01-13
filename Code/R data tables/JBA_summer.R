
library(tidyverse)
library(readxl)
library(here)
here::i_am(path = "./Code/R data tables/JBA_summer.R")


# seasons 
# JBA/AFBC Spring sampling (water data): March 30  - April 28 2021: 
# JBA/AFBC Summer sampling (water data): June 14 - July 16 2021
# JBA/AFBC Fall/Winter sampling (water data): November 19 2020 - January 14 2021

# WG/NRP Fall sampling (water data): October 9 - November 10 2020; second to last day 2020-11-05
# WG/NRP Summer sampling (water data): June 15 - August 19 (some sparse sampling) 2021
# WG/NRP Spring sampling (water data): February 9 - March 23 2021

# sample dates
# WG/NRP   as.Date("2020-11-10") (fall/winter)
            # as.Date("2021-04-28") (spring)
            # as.Date("2021-08-19") (summer) # <<< this file, median values
 
# JBA/AFBC as.Date("2021-04-28") (spring)
         # as.Date("2021-07-16") (summer)

# pull in relevant data:
d.fish<-read_xlsx(here("Data", "FINAL", "JBA_Biota_Data.xlsx"), sheet = "Data") #fish
d.water<-read_xlsx(here("Data", "FINAL", "JBA_Water_Data.xlsx"), sheet = "Data") #water
d.sed<-read_xlsx(here("Data", "FINAL", "JBA_Sediment_Data.xlsx"), sheet = "Data") #sediment

names(d.water)<-gsub(" ", "", names(d.water), fixed = TRUE)
names(d.sed)<-gsub(" ", "", names(d.sed), fixed = TRUE)
names(d.fish)<-gsub(" ", "", names(d.fish), fixed = TRUE)


# filter out data needed and format. 
d.fish<-d.fish %>%
  mutate_at(vars(Species, Location, Tissue), factor) %>% 
  select(Species,CommonName, Location, SampleDate, Tissue, Weight, 
         PFBA, PFPeA, PFHxA, PFHpA, PFOA, PFOS,
         PFNA, PFNS, PFDS,
         PFDA, PFDoS, PFDoA, PFDS,
         PFHpA, PFHxA, PFHxDA, PFHxS,
         PFTrDA, PFTeDA, PFPeS) %>% 
  mutate(Seasonality = if_else(SampleDate > as.Date("2021-06-01"), "summer", 
                            if_else(SampleDate < as.Date("2020-11-20"), "fall/winter",
                                    "spring"))) %>% 
  filter(Seasonality == "summer" &
           c(Tissue == "Whole Body" | Tissue == "Muscle")) %>% 
  mutate(PFUA = 0) %>% 
  mutate(Sp = if_else(CommonName == "Banded Killifish", "Kil",
                if_else(CommonName == "Creek Chubsucker", "Chu", 
                   if_else(CommonName == "Dace. Sp", "Dac",
                      if_else(CommonName == "Darter sp.", "Dar",
                        if_else(CommonName == "Eastern Mudminnow", "Min",
                          if_else(CommonName == "Margined Madtom", "Mad",
                            if_else(CommonName == "Pumpkinseed", "Pum",
                              if_else(CommonName == "Swallowtail Shiner", "Swa",
                                if_else(CommonName == "Fallfish","Fal",
                                  if_else(CommonName == "Large Mouth Bass", "Bas",
                                    if_else(CommonName == "Prey", "Pry", NA))))))))))))
# Fish masses
d.fish.w <- d.fish %>% 
  group_by(Sp, CommonName) %>% 
  summarise(n = n(), 
            median_mass = median(Weight, na.rm = TRUE), 
            min_mass = min(Weight, na.rm = TRUE), 
            max_mass = max(Weight, na.rm = TRUE))


# fillet-to-whole body conversion 
d.fish[d.fish$Tissue == "Muscle", c(6:22)]<-d.fish[d.fish$Tissue == "Muscle", c(6:22)]*2.5
d.fish<-d.fish %>% 
  dplyr::group_by(Sp) %>% 
  summarize(PFHxS.m = median(PFHxS), PFHxS.min = min(PFHxS), PFHxS.max = max(PFHxS),
            PFOS.m = median(PFOS), PFOS.min = min(PFOS), PFOS.max = max(PFOS),
            PFOA.m = median(PFOA), PFOA.min = min(PFOA), PFOA.max = max(PFOA),
            PFNA.m = median(PFNA), PFNA.min = min(PFNA), PFNA.max = max(PFNA),
            PFDA.m = median(PFDA, na.rm = TRUE),
            PFDA.min = min(PFDA, na.rm = TRUE),
            PFDA.max = max(PFDA, na.rm = TRUE),
            PFUA.m = median(PFUA), PFUA.min = min(PFUA), PFUA.max = max(PFUA))

d.water<-d.water %>%
  mutate_at(vars(SampleID, Location, Seasonality), factor) %>% 
  mutate(Seasonality = if_else(SampleDate > as.Date("2021-06-01"), "summer", 
                            if_else(SampleDate < as.Date("2020-11-20"), "fall/winter",
                                    "spring"))) %>% 
  select(SampleDate, Location, Seasonality, 
         Temp, pH, Cond, DO,
         PFBA, PFPeA, PFHxA, PFHpA, PFOA, PFOS,
         PFNA, PFNS, PFDS,
         PFDA, PFDoS, PFDoA, PFDS,
         PFHpA, PFHxA, PFHxDA, PFHxS,
         PFTrDA, PFTeDA, PFPeS) %>% 
  filter(Seasonality == "summer") %>% 
  mutate(PFUA = 0) %>% 
  summarize(PFHxS.m = median(PFHxS), PFHxS.min = min(PFHxS), PFHxS.max = max(PFHxS),
            PFOS.m = median(PFOS), PFOS.min = min(PFOS), PFOS.max = max(PFOS),
            PFOA.m = median(PFOA), PFOA.min = min(PFOA), PFOA.max = max(PFOA),
            PFNA.m = median(PFNA), PFNA.min = min(PFNA), PFNA.max = max(PFNA),
            PFDA.m = median(PFDA, na.rm = TRUE),
            PFDA.min = min(PFDA, na.rm = TRUE),
            PFDA.max = max(PFDA, na.rm = TRUE),
            PFUA.m = median(PFUA), PFUA.min = min(PFUA), PFUA.max = max(PFUA),
            temp.m = median(Temp), temp.min = min(Temp), temp.max = max(Temp),
            DO.m = median(DO, na.rm = TRUE),
            DO.min = min(DO, na.rm = TRUE),
            DO.max = max(DO, na.rm = TRUE))

d.sed<-d.sed%>%
  mutate_at(vars(SampleID, Location, Seasonality), factor) %>% 
  mutate(Seasonality = if_else(SampleDate > as.Date("2021-06-01"), "summer", 
                            if_else(SampleDate < as.Date("2020-11-20"), "fall/winter",
                                    "spring"))) %>% 
  select(SampleDate, Location, Seasonality, 
         PFBA, PFPeA, PFHxA, PFHpA, PFOA, PFOS,
         PFNA, PFNS, PFDS,
         PFDA, PFDoS, PFDoA, PFDS,
         PFHpA, PFHxA, PFHxS, 
         PFTrDA, PFTeDA, PFPeS) %>% 
  filter(Seasonality == "summer") %>% 
  mutate(PFUA = 0) %>% 
  summarize(PFHxS.m = median(PFHxS), PFHxS.min = min(PFHxS), PFHxS.max = max(PFHxS),
            PFOS.m = median(PFOS), PFOS.min = min(PFOS), PFOS.max = max(PFOS),
            PFOA.m = median(PFOA), PFOA.min = min(PFOA), PFOA.max = max(PFOA),
            PFNA.m = median(PFNA), PFNA.min = min(PFNA), PFNA.max = max(PFNA),
            PFDA.m = median(PFDA, na.rm = TRUE),
            PFDA.min = min(PFDA, na.rm = TRUE),
            PFDA.max = max(PFDA, na.rm = TRUE),
            PFUA.m = median(PFUA), PFUA.min = min(PFUA), PFUA.max = max(PFUA))




# order for PFAA is c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")
inputFiles_list<-create_data_tables(
  species = c("Phy",	"Bas",	"Swa",	"Dac",	"Min",	"Fal",	"Mad",	"Dar",	"Pum"), 
  group_species = c("plant", "fish", "fish", "fish", "fish", "fish", "fish", "fish", "fish"), 
  WB_kg = unlist(c(NA, 
            d.fish.w[d.fish.w$Sp == "Bas", "median_mass"]/1000,
            d.fish.w[d.fish.w$Sp == "Swa", "median_mass"]/1000, # in kg
            d.fish.w[d.fish.w$Sp == "Dac", "median_mass"]/1000,
            d.fish.w[d.fish.w$Sp == "Min", "median_mass"]/1000,
            d.fish.w[d.fish.w$Sp == "Fal", "median_mass"]/1000,
            d.fish.w[d.fish.w$Sp == "Mad", "median_mass"]/1000,
            d.fish.w[d.fish.w$Sp == "Dar", "median_mass"]/1000,
            d.fish.w[d.fish.w$Sp == "Pum", "median_mass"]/1000)), 
  m_O = c(1, 1, 1, 1, 1, 1, 1, 1, 1),
  GRF = c(0.8, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150),
  P_B = c(0.5, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15), 
  diet = list(
    Bas	= unlist(c(d.fish[d.fish$Sp == "Bas", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Bas", "PFOS.m"],
            d.fish[d.fish$Sp == "Bas", "PFOA.m"],
            d.fish[d.fish$Sp == "Bas", "PFNA.m"],
            d.fish[d.fish$Sp == "Bas", "PFDA.m"],
            d.fish[d.fish$Sp == "Bas", "PFUA.m"])), 
    Swa	= unlist(c(d.fish[d.fish$Sp == "Swa", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Swa", "PFOS.m"],
            d.fish[d.fish$Sp == "Swa", "PFOA.m"],
            d.fish[d.fish$Sp == "Swa", "PFNA.m"],
            d.fish[d.fish$Sp == "Swa", "PFDA.m"],
            d.fish[d.fish$Sp == "Swa", "PFUA.m"])),
    Dac	= unlist(c(d.fish[d.fish$Sp == "Dac", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Dac", "PFOS.m"],
            d.fish[d.fish$Sp == "Dac", "PFOA.m"],
            d.fish[d.fish$Sp == "Dac", "PFNA.m"],
            d.fish[d.fish$Sp == "Dac", "PFDA.m"],
            d.fish[d.fish$Sp == "Dac", "PFUA.m"])),
    Min	= unlist(c(d.fish[d.fish$Sp == "Min", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Min", "PFOS.m"],
            d.fish[d.fish$Sp == "Min", "PFOA.m"],
            d.fish[d.fish$Sp == "Min", "PFNA.m"],
            d.fish[d.fish$Sp == "Min", "PFDA.m"],
            d.fish[d.fish$Sp == "Min", "PFUA.m"])), #
    Fal	= unlist(c(d.fish[d.fish$Sp == "Fal", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Fal", "PFOS.m"],
            d.fish[d.fish$Sp == "Fal", "PFOA.m"],
            d.fish[d.fish$Sp == "Fal", "PFNA.m"],
            d.fish[d.fish$Sp == "Fal", "PFDA.m"],
            d.fish[d.fish$Sp == "Fal", "PFUA.m"])),
    Mad	= unlist(c(d.fish[d.fish$Sp == "Mad", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Mad", "PFOS.m"],
            d.fish[d.fish$Sp == "Mad", "PFOA.m"],
            d.fish[d.fish$Sp == "Mad", "PFNA.m"],
            d.fish[d.fish$Sp == "Mad", "PFDA.m"],
            d.fish[d.fish$Sp == "Mad", "PFUA.m"])),
    Dar	= unlist(c(d.fish[d.fish$Sp == "Dar", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Dar", "PFOS.m"],
            d.fish[d.fish$Sp == "Dar", "PFOA.m"],
            d.fish[d.fish$Sp == "Dar", "PFNA.m"],
            d.fish[d.fish$Sp == "Dar", "PFDA.m"],
            d.fish[d.fish$Sp == "Dar", "PFUA.m"])), # same
    Pum	= unlist(c(d.fish[d.fish$Sp == "Pum", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Pum", "PFOS.m"],
            d.fish[d.fish$Sp == "Pum", "PFOA.m"],
            d.fish[d.fish$Sp == "Pum", "PFNA.m"],
            d.fish[d.fish$Sp == "Pum", "PFDA.m"],
            d.fish[d.fish$Sp == "Pum", "PFUA.m"]))
      ),
  min_diet = list(
    Bas	= unlist(c(d.fish[d.fish$Sp == "Bas", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Bas", "PFOS.min"],
            d.fish[d.fish$Sp == "Bas", "PFOA.min"],
            d.fish[d.fish$Sp == "Bas", "PFNA.min"],
            d.fish[d.fish$Sp == "Bas", "PFDA.min"],
            d.fish[d.fish$Sp == "Bas", "PFUA.min"])), 
    Swa	= unlist(c(d.fish[d.fish$Sp == "Swa", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Swa", "PFOS.min"],
            d.fish[d.fish$Sp == "Swa", "PFOA.min"],
            d.fish[d.fish$Sp == "Swa", "PFNA.min"],
            d.fish[d.fish$Sp == "Swa", "PFDA.min"],
            d.fish[d.fish$Sp == "Swa", "PFUA.min"])),
    Dac	= unlist(c(d.fish[d.fish$Sp == "Dac", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Dac", "PFOS.min"],
            d.fish[d.fish$Sp == "Dac", "PFOA.min"],
            d.fish[d.fish$Sp == "Dac", "PFNA.min"],
            d.fish[d.fish$Sp == "Dac", "PFDA.min"],
            d.fish[d.fish$Sp == "Dac", "PFUA.min"])),
    Min	= unlist(c(d.fish[d.fish$Sp == "Min", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Min", "PFOS.min"],
            d.fish[d.fish$Sp == "Min", "PFOA.min"],
            d.fish[d.fish$Sp == "Min", "PFNA.min"],
            d.fish[d.fish$Sp == "Min", "PFDA.min"],
            d.fish[d.fish$Sp == "Min", "PFUA.min"])), #
    Fal	= unlist(c(d.fish[d.fish$Sp == "Fal", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Fal", "PFOS.min"],
            d.fish[d.fish$Sp == "Fal", "PFOA.min"],
            d.fish[d.fish$Sp == "Fal", "PFNA.min"],
            d.fish[d.fish$Sp == "Fal", "PFDA.min"],
            d.fish[d.fish$Sp == "Fal", "PFUA.min"])),
    Mad	= unlist(c(d.fish[d.fish$Sp == "Mad", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Mad", "PFOS.min"],
            d.fish[d.fish$Sp == "Mad", "PFOA.min"],
            d.fish[d.fish$Sp == "Mad", "PFNA.min"],
            d.fish[d.fish$Sp == "Mad", "PFDA.min"],
            d.fish[d.fish$Sp == "Mad", "PFUA.min"])),
    Dar	= unlist(c(d.fish[d.fish$Sp == "Dar", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Dar", "PFOS.min"],
            d.fish[d.fish$Sp == "Dar", "PFOA.min"],
            d.fish[d.fish$Sp == "Dar", "PFNA.min"],
            d.fish[d.fish$Sp == "Dar", "PFDA.min"],
            d.fish[d.fish$Sp == "Dar", "PFUA.min"])), # same
    Pum	= unlist(c(d.fish[d.fish$Sp == "Pum", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Pum", "PFOS.min"],
            d.fish[d.fish$Sp == "Pum", "PFOA.min"],
            d.fish[d.fish$Sp == "Pum", "PFNA.min"],
            d.fish[d.fish$Sp == "Pum", "PFDA.min"],
            d.fish[d.fish$Sp == "Pum", "PFUA.min"]))),
  max_diet = list(
    Bas	= unlist(c(d.fish[d.fish$Sp == "Bas", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Bas", "PFOS.max"],
            d.fish[d.fish$Sp == "Bas", "PFOA.max"],
            d.fish[d.fish$Sp == "Bas", "PFNA.max"],
            d.fish[d.fish$Sp == "Bas", "PFDA.max"],
            d.fish[d.fish$Sp == "Bas", "PFUA.max"])), 
    Swa	= unlist(c(d.fish[d.fish$Sp == "Swa", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Swa", "PFOS.max"],
            d.fish[d.fish$Sp == "Swa", "PFOA.max"],
            d.fish[d.fish$Sp == "Swa", "PFNA.max"],
            d.fish[d.fish$Sp == "Swa", "PFDA.max"],
            d.fish[d.fish$Sp == "Swa", "PFUA.max"])),
    Dac	= unlist(c(d.fish[d.fish$Sp == "Dac", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Dac", "PFOS.max"],
            d.fish[d.fish$Sp == "Dac", "PFOA.max"],
            d.fish[d.fish$Sp == "Dac", "PFNA.max"],
            d.fish[d.fish$Sp == "Dac", "PFDA.max"],
            d.fish[d.fish$Sp == "Dac", "PFUA.max"])),
    Min	= unlist(c(d.fish[d.fish$Sp == "Min", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Min", "PFOS.max"],
            d.fish[d.fish$Sp == "Min", "PFOA.max"],
            d.fish[d.fish$Sp == "Min", "PFNA.max"],
            d.fish[d.fish$Sp == "Min", "PFDA.max"],
            d.fish[d.fish$Sp == "Min", "PFUA.max"])), #
    Fal	= unlist(c(d.fish[d.fish$Sp == "Fal", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Fal", "PFOS.max"],
            d.fish[d.fish$Sp == "Fal", "PFOA.max"],
            d.fish[d.fish$Sp == "Fal", "PFNA.max"],
            d.fish[d.fish$Sp == "Fal", "PFDA.max"],
            d.fish[d.fish$Sp == "Fal", "PFUA.max"])),
    Mad	= unlist(c(d.fish[d.fish$Sp == "Mad", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Mad", "PFOS.max"],
            d.fish[d.fish$Sp == "Mad", "PFOA.max"],
            d.fish[d.fish$Sp == "Mad", "PFNA.max"],
            d.fish[d.fish$Sp == "Mad", "PFDA.max"],
            d.fish[d.fish$Sp == "Mad", "PFUA.max"])),
    Dar	= unlist(c(d.fish[d.fish$Sp == "Dar", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Dar", "PFOS.max"],
            d.fish[d.fish$Sp == "Dar", "PFOA.max"],
            d.fish[d.fish$Sp == "Dar", "PFNA.max"],
            d.fish[d.fish$Sp == "Dar", "PFDA.max"],
            d.fish[d.fish$Sp == "Dar", "PFUA.max"])), # same
    Pum	= unlist(c(d.fish[d.fish$Sp == "Pum", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Pum", "PFOS.max"],
            d.fish[d.fish$Sp == "Pum", "PFOA.max"],
            d.fish[d.fish$Sp == "Pum", "PFNA.max"],
            d.fish[d.fish$Sp == "Pum", "PFDA.max"],
            d.fish[d.fish$Sp == "Pum", "PFUA.max"]))),
  foodWeb = list(
Phy = c(0,0,0,0,0,0,0,0,0,0),
Bas = c(0.8,0.2,0,0,0,0,0,0,0,0),
Swa = c(0.3,0.6,0,0,0,0,0,0,0,0),
Dac = c(0.5,0.4,0,0,0,0,0,0,0,0),
Min = c(0.4,0.6,0,0,0,0,0,0,0,0),
Fal = c(0.2,0.3,0,0.1,0.1,0.1,0,0,0,0),
Mad = c(0.2,0.3,0,0.1,0.1,0.1,0,0,0,0),
Dar = c(0.2,0.3,0,0.1,0,0,0,0,0,0),
Pum = c(0.2,0.2,0.1,0,0.1,0.1,0,0,0,0)
      ),
  C_WTO_ng_mL = unlist(c(d.water[, "PFHxS.m"]/1000, # in ng/mL
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
  C_s_ng_g = unlist(c(d.sed[, "PFHxS.m"], 
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
  C_OX = unlist(c(d.water[, "DO.m"])), 
  T = unlist(c(d.water[, "temp.m"])),
  # log_Koc = c("calibr", 2.85	, 2.9, 2.203283149, 2.932949081, 5,	4.9900)
  log_Koc = c("Brown etal", 2.3317871	, 2.891004, 2.3032831,	3.0329491,	6,	5) # c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")
  
# ,
  # log_Dmw = c("calibr-cosmotherm", 3.37, 4.61, 3.47, 4.10, 4.69, 5.34)


) # order for PFAA)

