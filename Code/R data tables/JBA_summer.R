
library(tidyverse)
library(readxl)
library(here)
here::i_am(path = "./Code/R data tables/JBA_summer.R")


season = "Summer"
# pull in relevant data:
d.fish<-read.csv(here("Data", "Original_Brown_etal", "JBA", "JBA_Biota_Data.csv")) #fish
d.fish.w<-read_xlsx(path = here("Data", "Original_Brown_etal", "JBA", "JBA_Biota_fishmass_kk.xlsx"),
                    sheet = "data")
d.water<-read.csv(here("Data", "Original_Brown_etal", "JBA", "JBA_Water_Data.csv")) #water
d.sed<-read.csv(here("Data", "Original_Brown_etal", "JBA", "JBA_Sediment_Data.csv")) #sediment

# get sample sizes for everything first
n.fish<-d.fish %>% 
  dplyr::group_by(Seasonality, Species) %>% 
  filter(Tissue == "Muscle" | Tissue == "Whole Body") %>% 
  # mutate(sampleID2 = substr(SampleID, start = 1, stop = 6)) %>% 
  summarize(n = n())
  
n.sed<-d.sed %>% 
  mutate(sampleID2 = substr(Sample.ID, start = 1, stop = 6)) %>% 
  summarize(n = n_distinct(sampleID2))

n.water<-d.water %>% # temperature and PFAA
  mutate(sampleID2 = substr(SampleID, start = 1, stop = 6)) %>% 
  summarize(n = n_distinct(sampleID2))

# filter out data needed and format. 
d.fish<-d.fish %>%
  mutate_at(vars(Species, Location, Seasonality, Tissue), factor) %>% 
  select(Species, Location, Seasonality, Tissue, Length, 
         PFBA, PFPeA, PFHxA, PFHpA, PFOA, PFOS,
         PFNA, PFNS, PFDS,
         PFDA, PFDoS, PFDoA, PFDS,
         PFHpA, PFHxA, PFHxDA, PFHxS,
         PFTrDA, PFTeDA, PFPeS, PFUdA) %>% 
  filter(Seasonality == season &
           c(Tissue == "Whole Body" | Tissue == "Muscle")) %>% 
  mutate(PFUA = 0) %>% 
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
                                    if_else(Species == "Prey", "Pry", NA))))))))))))

# fillet-to-whole body conversion 
d.fish[d.fish$Tissue == "Muscle", c(6:22)]<-d.fish[d.fish$Tissue == "Muscle", c(6:22)]*2.5
d.fish<-d.fish %>% 
  dplyr::group_by(Sp) %>% 
  summarize(PFHxS.m = median(PFHxS), PFHxS.min = min(PFHxS), PFHxS.max = max(PFHxS),
            PFOS.m = median(PFOS), PFOS.min = min(PFOS), PFOS.max = max(PFOS),
            PFOA.m = median(PFOA), PFOA.min = min(PFOA), PFOA.max = max(PFOA),
            PFNA.m = median(PFNA), PFNA.min = min(PFNA), PFNA.max = max(PFNA),
            PFDA.m = median(PFDA), PFDA.min = min(PFDA), PFDA.max = max(PFDA),
            PFUA.m = median(PFUA), PFUA.min = min(PFUA), PFUA.max = max(PFUA))

d.water<-d.water %>%
  mutate_at(vars(SampleID, Location, Seasonality), factor) %>% 
  select(SampleDate, Location, Seasonality, 
         Temp, pH, Cond, DO,
         PFBA, PFPeA, PFHxA, PFHpA, PFOA, PFOS,
         PFNA, PFNS, PFDS,
         PFDA, PFDoS, PFDoA, PFDS,
         PFHpA, PFHxA, PFHxDA, PFHxS,
         PFTrDA, PFTeDA, PFPeS) %>% 
  filter(Seasonality == season) %>% 
  mutate(PFUA = 0) %>% 
  dplyr::group_by(Seasonality) %>% 
  summarize(PFHxS.m = median(PFHxS), PFHxS.min = min(PFHxS), PFHxS.max = max(PFHxS),
            PFOS.m = median(PFOS), PFOS.min = min(PFOS), PFOS.max = max(PFOS),
            PFOA.m = median(PFOA), PFOA.min = min(PFOA), PFOA.max = max(PFOA),
            PFNA.m = median(PFNA), PFNA.min = min(PFNA), PFNA.max = max(PFNA),
            PFDA.m = median(PFDA), PFDA.min = min(PFDA), PFDA.max = max(PFDA),
            PFUA.m = median(PFUA), PFUA.min = min(PFUA), PFUA.max = max(PFUA),
            temp.m = median(Temp), temp.min = min(Temp), temp.max = max(Temp),
            DO.m = median(DO), DO.min = min(DO), DO.max = max(DO))

d.sed<-d.sed%>%
  mutate_at(vars(Sample.ID, Location, Seasonality), factor) %>% 
  select(SampleDate, Location, Seasonality, 
         PFBA, PFPeA, PFHxA, PFHpA, PFOA, PFOS,
         PFNA, PFNS, PFDS,
         PFDA, PFDoS, PFDoA, PFDS,
         PFHpA, PFHxA, PFHxS, 
         PFTrDA, PFTeDA, PFPeS) %>% 
  filter(Seasonality == season) %>% 
  mutate(PFUA = 0) %>% 
  dplyr::group_by(Seasonality) %>% 
  summarize(PFHxS.m = median(PFHxS), PFHxS.min = min(PFHxS), PFHxS.max = max(PFHxS),
            PFOS.m = median(PFOS), PFOS.min = min(PFOS), PFOS.max = max(PFOS),
            PFOA.m = median(PFOA), PFOA.min = min(PFOA), PFOA.max = max(PFOA),
            PFNA.m = median(PFNA), PFNA.min = min(PFNA), PFNA.max = max(PFNA),
            PFDA.m = median(PFDA), PFDA.min = min(PFDA), PFDA.max = max(PFDA),
            PFUA.m = median(PFUA), PFUA.min = min(PFUA), PFUA.max = max(PFUA))


# order for PFAA is c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")
inputFiles_list<-create_data_tables(
  species = c("Phy",	"Dac",	"Dar",	"Min",	"Fal",	"Bas",	"Mad",	"Pum",	"Swa"), 
  group_species = c("plant", "fish", "fish", "fish", "fish", "fish", "fish", "fish", "fish"), 
  WB_kg = unlist(c(NA, 
            d.fish.w[d.fish.w$Sp == "Dac", "mean_mass"]/1000, # in kg
            d.fish.w[d.fish.w$Sp == "Dar", "mean_mass"]/1000,
            d.fish.w[d.fish.w$Sp == "Min", "mean_mass"]/1000,
            d.fish.w[d.fish.w$Sp == "Fal", "mean_mass"]/1000,
            d.fish.w[d.fish.w$Sp == "Bas", "mean_mass"]/1000,
            d.fish.w[d.fish.w$Sp == "Mad", "mean_mass"]/1000,
            d.fish.w[d.fish.w$Sp == "Pum", "mean_mass"]/1000,
            d.fish.w[d.fish.w$Sp == "Swa", "mean_mass"]/1000)), 
  m_O = c(1, 1, 1, 1, 1, 1, 1, 1, 1),
  GRF = c(0.8, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150),
  P_B = c(0.5, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15), 
  diet = list(
    Dac	= unlist(c(d.fish[d.fish$Sp == "Dac", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Dac", "PFOS.m"],
            d.fish[d.fish$Sp == "Dac", "PFOA.m"],
            d.fish[d.fish$Sp == "Dac", "PFNA.m"],
            d.fish[d.fish$Sp == "Dac", "PFDA.m"],
            d.fish[d.fish$Sp == "Dac", "PFUA.m"])),
    Dar	= unlist(c(d.fish[d.fish$Sp == "Dar", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Dar", "PFOS.m"],
            d.fish[d.fish$Sp == "Dar", "PFOA.m"],
            d.fish[d.fish$Sp == "Dar", "PFNA.m"],
            d.fish[d.fish$Sp == "Dar", "PFDA.m"],
            d.fish[d.fish$Sp == "Dar", "PFUA.m"])), # same
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
    Bas	= unlist(c(d.fish[d.fish$Sp == "Bas", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Bas", "PFOS.m"],
            d.fish[d.fish$Sp == "Bas", "PFOA.m"],
            d.fish[d.fish$Sp == "Bas", "PFNA.m"],
            d.fish[d.fish$Sp == "Bas", "PFDA.m"],
            d.fish[d.fish$Sp == "Bas", "PFUA.m"])), 
    Mad	= unlist(c(d.fish[d.fish$Sp == "Mad", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Mad", "PFOS.m"],
            d.fish[d.fish$Sp == "Mad", "PFOA.m"],
            d.fish[d.fish$Sp == "Mad", "PFNA.m"],
            d.fish[d.fish$Sp == "Mad", "PFDA.m"],
            d.fish[d.fish$Sp == "Mad", "PFUA.m"])),
    Pum	= unlist(c(d.fish[d.fish$Sp == "Pum", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Pum", "PFOS.m"],
            d.fish[d.fish$Sp == "Pum", "PFOA.m"],
            d.fish[d.fish$Sp == "Pum", "PFNA.m"],
            d.fish[d.fish$Sp == "Pum", "PFDA.m"],
            d.fish[d.fish$Sp == "Pum", "PFUA.m"])),
    Swa	= unlist(c(d.fish[d.fish$Sp == "Swa", "PFHxS.m"], # in ng/g
            d.fish[d.fish$Sp == "Swa", "PFOS.m"],
            d.fish[d.fish$Sp == "Swa", "PFOA.m"],
            d.fish[d.fish$Sp == "Swa", "PFNA.m"],
            d.fish[d.fish$Sp == "Swa", "PFDA.m"],
            d.fish[d.fish$Sp == "Swa", "PFUA.m"])
      )),
  min_diet = list(
    Dac	= unlist(c(d.fish[d.fish$Sp == "Dac", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Dac", "PFOS.min"],
            d.fish[d.fish$Sp == "Dac", "PFOA.min"],
            d.fish[d.fish$Sp == "Dac", "PFNA.min"],
            d.fish[d.fish$Sp == "Dac", "PFDA.min"],
            d.fish[d.fish$Sp == "Dac", "PFUA.min"])),
    Dar	= unlist(c(d.fish[d.fish$Sp == "Dar", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Dar", "PFOS.min"],
            d.fish[d.fish$Sp == "Dar", "PFOA.min"],
            d.fish[d.fish$Sp == "Dar", "PFNA.min"],
            d.fish[d.fish$Sp == "Dar", "PFDA.min"],
            d.fish[d.fish$Sp == "Dar", "PFUA.min"])), # same
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
    Bas	= unlist(c(d.fish[d.fish$Sp == "Bas", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Bas", "PFOS.min"],
            d.fish[d.fish$Sp == "Bas", "PFOA.min"],
            d.fish[d.fish$Sp == "Bas", "PFNA.min"],
            d.fish[d.fish$Sp == "Bas", "PFDA.min"],
            d.fish[d.fish$Sp == "Bas", "PFUA.min"])), 
    Mad	= unlist(c(d.fish[d.fish$Sp == "Mad", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Mad", "PFOS.min"],
            d.fish[d.fish$Sp == "Mad", "PFOA.min"],
            d.fish[d.fish$Sp == "Mad", "PFNA.min"],
            d.fish[d.fish$Sp == "Mad", "PFDA.min"],
            d.fish[d.fish$Sp == "Mad", "PFUA.min"])),
    Pum	= unlist(c(d.fish[d.fish$Sp == "Pum", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Pum", "PFOS.min"],
            d.fish[d.fish$Sp == "Pum", "PFOA.min"],
            d.fish[d.fish$Sp == "Pum", "PFNA.min"],
            d.fish[d.fish$Sp == "Pum", "PFDA.min"],
            d.fish[d.fish$Sp == "Pum", "PFUA.min"])),
    Swa	= unlist(c(d.fish[d.fish$Sp == "Swa", "PFHxS.min"], # in ng/g
            d.fish[d.fish$Sp == "Swa", "PFOS.min"],
            d.fish[d.fish$Sp == "Swa", "PFOA.min"],
            d.fish[d.fish$Sp == "Swa", "PFNA.min"],
            d.fish[d.fish$Sp == "Swa", "PFDA.min"],
            d.fish[d.fish$Sp == "Swa", "PFUA.min"])
      )),
  max_diet = list(
    Dac	= unlist(c(d.fish[d.fish$Sp == "Dac", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Dac", "PFOS.max"],
            d.fish[d.fish$Sp == "Dac", "PFOA.max"],
            d.fish[d.fish$Sp == "Dac", "PFNA.max"],
            d.fish[d.fish$Sp == "Dac", "PFDA.max"],
            d.fish[d.fish$Sp == "Dac", "PFUA.max"])),
    Dar	= unlist(c(d.fish[d.fish$Sp == "Dar", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Dar", "PFOS.max"],
            d.fish[d.fish$Sp == "Dar", "PFOA.max"],
            d.fish[d.fish$Sp == "Dar", "PFNA.max"],
            d.fish[d.fish$Sp == "Dar", "PFDA.max"],
            d.fish[d.fish$Sp == "Dar", "PFUA.max"])), # same
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
    Bas	= unlist(c(d.fish[d.fish$Sp == "Bas", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Bas", "PFOS.max"],
            d.fish[d.fish$Sp == "Bas", "PFOA.max"],
            d.fish[d.fish$Sp == "Bas", "PFNA.max"],
            d.fish[d.fish$Sp == "Bas", "PFDA.max"],
            d.fish[d.fish$Sp == "Bas", "PFUA.max"])), 
    Mad	= unlist(c(d.fish[d.fish$Sp == "Mad", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Mad", "PFOS.max"],
            d.fish[d.fish$Sp == "Mad", "PFOA.max"],
            d.fish[d.fish$Sp == "Mad", "PFNA.max"],
            d.fish[d.fish$Sp == "Mad", "PFDA.max"],
            d.fish[d.fish$Sp == "Mad", "PFUA.max"])),
    Pum	= unlist(c(d.fish[d.fish$Sp == "Pum", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Pum", "PFOS.max"],
            d.fish[d.fish$Sp == "Pum", "PFOA.max"],
            d.fish[d.fish$Sp == "Pum", "PFNA.max"],
            d.fish[d.fish$Sp == "Pum", "PFDA.max"],
            d.fish[d.fish$Sp == "Pum", "PFUA.max"])),
    Swa	= unlist(c(d.fish[d.fish$Sp == "Swa", "PFHxS.max"], # in ng/g
            d.fish[d.fish$Sp == "Swa", "PFOS.max"],
            d.fish[d.fish$Sp == "Swa", "PFOA.max"],
            d.fish[d.fish$Sp == "Swa", "PFNA.max"],
            d.fish[d.fish$Sp == "Swa", "PFDA.max"],
            d.fish[d.fish$Sp == "Swa", "PFUA.min"])
      )),
  foodWeb = list(
    # Phy	= c(0,	  0,	  0,	  0,	  0,	  0,	0,	0,	0,	0),
    # Dac	= c(0.5,	0.5,	0,	  0,	  0,  	0,	0,	0,	0,	0),
    # Dar	= c(0.1,	0,	  0.1,	0,	  0.1,	0,	0,	0,	0,	0.1),
    # Min	= c(0.4,	0.5,	0,	  0,	  0,	  0,	0,	0,	0,	0.1),
    # Fal	= c(0,	  0,	  0.15,	0.15,	0.5,	0,	0,	0,	0,	0.1),
    # Bas	= c(0,	  0,	  0.1,	0.1,	0,    0,	0,	0,	0.5,	0.1),
    # Mad	= c(0.5,	0.5,	0,	  0,	  0,	  0,	0,	0,	0,	0),
    # Pum	= c(0.1,	0.1,	0.1,	0.1,	0,	  0,	0,	0,	0,	0.1),
    # Swa	= c(0.5,	0.5,	0,	  0,	  0,	  0,	0,	0,	0,	0)
Phy = c(1,0,0,0,0,0,0,0,0,0),
Dac = c(0.5,0.4,0,0,0,0,0,0,0,0),
Dar = c(0.5,0.6,0,0,0,0,0,0,0,0),
Min = c(0.4,0.6,0,0,0,0,0,0,0,0),
Fal = c(0.5,0.5,0,0,0,0,0,0,0,0),
Bas = c(0.2,0.4,0,0.1,0.1,0,0,0,0,0.2),
Mad = c(0.5,0.5,0,0,0,0,0,0,0,0),
Pum = c(0.3,0.5,0.1,0,0.05,0,0,0,0,0.05),
Swa = c(0.3,0.7,0,0,0,0,0,0,0,0)
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
  log_Koc = c("calibr", 2.85	, 2.9, 2.203283149,	6.983623,6,	4.9900)
  
# ,
  # log_Dmw = c("calibr-cosmotherm", 3.37, 4.61, 3.47, 4.10, 4.69, 5.34)


) # order for PFAA)

