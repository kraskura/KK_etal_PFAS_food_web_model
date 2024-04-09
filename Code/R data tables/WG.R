
library(tidyverse)
library(readxl)
library(here)
here::i_am(path = "./Code/R data tables/WG.R")

# pull in relevant data:
d.fish<-read.csv(here("Data", "Original_Brown_etal", "WG", "WG_Biota_Data.csv")) #fish
d.fish.w<-read_xlsx(path = here("Data", "Original_Brown_etal", "WG", "WG_fishMass_kk.xlsx"),
                    sheet = "data")
d.water<-read.csv(here("Data", "Original_Brown_etal", "WG", "WG_Water_Data.csv")) #water
d.sed<-read.csv(here("Data", "Original_Brown_etal", "WG", "WG_Sediment_Data.csv")) #sediment

# get sample sizes for everything first
n.fish<-d.fish %>% 
  mutate(sampleID2 = substr(SampleID, start = 1, stop = 6)) %>% 
  dplyr:::group_by(Species) %>% 
  summarize(n = n_distinct(sampleID2))
  
n.sed<-d.sed %>% 
  mutate(sampleID2 = substr(SampleID, start = 1, stop = 6)) %>% 
  summarize(n = n_distinct(sampleID2))

n.water<-d.water %>% # temperature and PFAA
  mutate(sampleID2 = substr(SampleID, start = 1, stop = 6)) %>% 
  summarize(n = n_distinct(sampleID2))

# filter out data needed and format. 
d.fish<-d.fish %>%
  mutate_at(vars(Species, Tissue), factor) %>% 
  select(Species, Tissue, Length, 
         PFBA, PFPeA, PFHxA, PFHpA, PFOA, PFOS,
         PFNA, PFHpA, PFHxA, PFHxS, PFPeS) %>% 
  filter(c(Tissue == "WholeBody" | Tissue == "Muscle")) %>% 
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

# fillet-to-whole body conversion 
d.fish[d.fish$Tissue == "Muscle", c(4:14)]<-d.fish[d.fish$Tissue == "Muscle", c(4:14)]*2.5
d.fish<-d.fish %>% 
  dplyr::group_by(Sp) %>% 
  summarize(PFHxS.m = median(PFHxS), PFHxS.min = min(PFHxS), PFHxS.max = max(PFHxS),
            PFOS.m = median(PFOS), PFOS.min = min(PFOS), PFOS.max = max(PFOS),
            PFOA.m = median(PFOA), PFOA.min = min(PFOA), PFOA.max = max(PFOA),
            PFNA.m = median(PFNA), PFNA.min = min(PFNA), PFNA.max = max(PFNA),
            PFDA.m = median(PFDA), PFDA.min = min(PFDA), PFDA.max = max(PFDA),
            PFUA.m = median(PFUA), PFUA.min = min(PFUA), PFUA.max = max(PFUA))

d.water<-d.water %>%
  mutate_at(vars(SampleID, Location), factor) %>% 
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
            temp.m = median(WaterTemp), temp.min = min(WaterTemp), temp.max = max(WaterTemp))

d.sed<-d.sed%>%
  mutate_at(vars(SampleID, Location), factor) %>% 
  select(SampleDate, Location,
         PFBA, PFPeA, PFHxA, PFHpA, PFOA, PFOS,
         PFNA, PFHpA, PFHxA, PFHxS, PFPeS) %>% 
  mutate(PFUA = 0, PFDA = 0) %>% 
  summarize(PFHxS.m = median(PFHxS), PFHxS.min = min(PFHxS), PFHxS.max = max(PFHxS),
            PFOS.m = median(PFOS), PFOS.min = min(PFOS), PFOS.max = max(PFOS),
            PFOA.m = median(PFOA), PFOA.min = min(PFOA), PFOA.max = max(PFOA),
            PFNA.m = median(PFNA), PFNA.min = min(PFNA), PFNA.max = max(PFNA),
            PFDA.m = median(PFDA), PFDA.min = min(PFDA), PFDA.max = max(PFDA),
            PFUA.m = median(PFUA), PFUA.min = min(PFUA), PFUA.max = max(PFUA))


# order for PFAA is c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")
inputFiles_list<-create_data_tables(
  species = c("Phy", "Pry", "Bgl", "Bas"), 
  group_species = c("plant", "fish", "fish", "fish"), 
  WB_kg = unlist(c(NA, 
            d.fish.w[d.fish.w$Sp == "Pry", "mean_mass"]/1000, # in kg
            d.fish.w[d.fish.w$Sp == "Bgl", "mean_mass"]/1000,
            d.fish.w[d.fish.w$Sp == "Bas", "mean_mass"]/1000)), 
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
           			
) # order for PFAA)


