

# Title: Data prep for monte carlo simulation 
# Purpose: extracts the min, max, mean values for fish, water, and sediment  PFAS 
# concentrations, fish body sizes, and water temperatures. 
# These values are used in mc_sim_data_table.R to randomly select the mentioned 
# values when constructing inputFiles_list for a model run. 
# 
# library(readxl)
# library(tidyverse)
# 
# metadata<-read_xlsx("Data/Species_specific_info.xlsx",sheet =1)

library(tidyverse)
library(readxl)
library(here)
# # pull in relevant data:
d.fish.jba<-read.csv(here("Data/Original_Brown_etal/JBA/JBA_Biota_Data.csv")) #fish
d.fish.w.jba<-read_xlsx(here("Data/Original_Brown_etal/JBA/JBA_Biota_fishmass_kk.xlsx"), sheet = "data")
d.water.jba<-read.csv(here("Data/Original_Brown_etal/JBA/JBA_Water_Data.csv")) #water
d.sed.jba<-read.csv(here("Data/Original_Brown_etal/JBA/JBA_Sediment_Data.csv")) #sediment

d.fish.wg<-read.csv(here("Data/Original_Brown_etal/WG/WG_Biota_Data.csv")) #fish
d.fish.w.wg<-read_xlsx(here("Data/Original_Brown_etal/WG/WG_fishMass_kk.xlsx"), sheet = "data")
d.water.wg<-read.csv(here("Data/Original_Brown_etal/WG/WG_Water_Data.csv")) #water
d.sed.wg<-read.csv(here("Data/Original_Brown_etal/WG/WG_Sediment_Data.csv")) #sediment

d.fish.wg$Seasonality<-"Fall"
d.fish.w.wg$Seasonality<-"Fall"
d.water.wg$Seasonality<-"Fall"
d.sed.wg$Seasonality<-"Fall"

d.fish.wg <- d.fish.wg %>% 
    select(Species, Seasonality, Tissue, Length, 
         PFOA, PFOS,
         PFNA, PFHxS) 
d.fish.jba <- d.fish.jba %>% 
    select(Species, Seasonality, Tissue, Length, 
         PFOA, PFOS,
         PFNA, PFHxS)  

# filter out data needed and format. 
d.fish<-rbind(d.fish.wg, d.fish.jba) %>%
  mutate_at(vars(Species,  Seasonality, Tissue), factor) %>% 
  filter(c(Tissue == "WholeBody" | Tissue == "Whole Body" | Tissue == "Muscle")) %>% 
  mutate(PFUA = 0, PFDA = 0) %>% 
  mutate(loc = if_else(Seasonality == "Fall", "WG",
                if_else(Seasonality == "Spring", "JBA", 
                 if_else(Seasonality == "Summer", "JBA", NA)))) %>% 
  mutate(SppAlias = if_else(Species == "Banded Killifish", "Kil",
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
d.fish[d.fish$Tissue == "Muscle", c(5:8)]<-d.fish[d.fish$Tissue == "Muscle", c(5:8)]*2.5

# sample size of each measurement 
n_fish <- d.fish %>% 
  group_by(loc, Seasonality, SppAlias) %>% 
  dplyr::summarize(n = n(),
                   min_size = min(Length),
                   max_length = max(Length)) %>%
  as.data.frame()

# water dataset 
d.water.wg <- d.water.wg %>% 
  select(Seasonality, WaterTemp,
         PFOA, PFOS,
         PFNA, PFHxS) 
d.water.jba <- d.water.jba %>% 
  select(Seasonality, Temp,
         PFOA, PFOS,
         PFNA, PFHxS) %>% 
  dplyr:::rename(WaterTemp = Temp)

d.water<-rbind(d.water.wg, d.water.jba) %>%
  mutate_at(vars(Seasonality), factor) %>% 
  mutate(loc = if_else(Seasonality == "Fall", "WG",
                if_else(Seasonality == "Spring", "JBA", 
                 if_else(Seasonality == "Summer", "JBA", NA)))) %>%
  mutate(PFUA = 0, PFDA = 0)%>% 
  dplyr:::rename(Temperature = WaterTemp)

# sediment dataset 
d.sed.wg <- d.sed.wg %>% 
  select(Seasonality, 
         PFOA, PFOS,
         PFNA, PFHxS) 
d.sed.jba <- d.sed.jba %>% 
  select(Seasonality, 
         PFOA, PFOS,
         PFNA, PFHxS) 

d.sed <- rbind(d.sed.jba, d.sed.wg) %>%
  mutate_at(vars(Seasonality), factor) %>% 
  mutate(loc = if_else(Seasonality == "Fall", "WG",
                if_else(Seasonality == "Spring", "JBA", 
                 if_else(Seasonality == "Summer", "JBA", NA)))) %>%
  mutate(PFUA = 0, PFDA = 0) 
  
# means datasets ********************
# carcass only 
dfish.MEAN<-d.fish %>% 
  dplyr:::group_by(loc, SppAlias, Seasonality) %>% 
  summarize(PFHxS = median(PFHxS),
            PFOS =  median(PFOS), 
            PFOA =  median(PFOA), 
            PFNA =  median(PFNA),
            PFDA =  median(PFDA), 
            PFUA =  median(PFUA),.groups = 'keep') 
dfish.MIN<-d.fish %>% 
  dplyr:::group_by( loc, SppAlias, Seasonality) %>% 
  summarize(
            PFHxS = min(PFHxS),
            PFOS = min(PFOS), 
            PFOA = min(PFOA), 
            PFNA = min(PFNA),
            PFDA = min(PFDA), 
            PFUA = min(PFUA),.groups = 'keep') 
dfish.MAX<-d.fish %>% 
  dplyr:::group_by( loc, SppAlias, Seasonality) %>% 
  summarize(
            PFHxS = max(PFHxS),
            PFOS = max(PFOS), 
            PFOA = max(PFOA), 
            PFNA = max(PFNA),
            PFDA = max(PFDA), 
            PFUA = max(PFUA),.groups = 'keep') 

dfish.CV<-d.fish %>% 
  dplyr:::group_by( loc, SppAlias, Seasonality) %>% 
  summarize(
            PFHxS = sd(PFHxS)/mean(PFHxS),
            PFOS = sd(PFOS)/mean(PFOS), 
            PFOA = sd(PFOA)/mean(PFOA), 
            PFNA = sd(PFNA)/mean(PFNA),
            PFDA = sd(PFDA)/mean(PFDA), 
            PFUA = sd(PFUA)/mean(PFUA),.groups = 'keep') %>% 
  pivot_longer(cols = c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA"), 
               names_to = "PFAA", values_to = "fish_cv")

# water samples
dwater.MEAN<-d.water %>% 
  dplyr:::group_by(loc, Seasonality) %>% 
  summarize(PFHxS = median(PFHxS),
            PFOS = median(PFOS), 
            PFOA = median(PFOA), 
            PFNA = median(PFNA),
            PFDA = median(PFDA), 
            PFUA = median(PFUA), 
            Temp = median(Temperature),.groups = 'keep') 

dwater.MIN<-d.water  %>% 
  dplyr:::group_by(loc, Seasonality) %>%
  summarize(PFHxS = min(PFHxS),
            PFOS = min(PFOS), 
            PFOA = min(PFOA), 
            PFNA = min(PFNA),
            PFDA = min(PFDA), 
            PFUA = min(PFUA),
            Temp = min(Temperature),.groups = 'keep') 

dwater.MAX<-d.water  %>% 
  dplyr:::group_by(loc, Seasonality) %>% 
  summarize(
            PFHxS = max(PFHxS),
            PFOS = max(PFOS), 
            PFOA = max(PFOA), 
            PFNA = max(PFNA),
            PFDA = max(PFDA), 
            PFUA = max(PFUA), 
            Temp = max(Temperature),.groups = 'keep') 

# sediment samples
dsed.MEAN<-d.sed %>% 
  dplyr:::group_by(loc, Seasonality) %>% 
  summarize(PFHxS = median(PFHxS),
            PFOS = median(PFOS), 
            PFOA = median(PFOA), 
            PFNA = median(PFNA),
            PFDA = median(PFDA), 
            PFUA = median(PFUA), .groups = 'keep') 

dsed.MIN<-d.sed  %>% 
  dplyr:::group_by(loc, Seasonality) %>%
  summarize(PFHxS = min(PFHxS),
            PFOS = min(PFOS), 
            PFOA = min(PFOA), 
            PFNA = min(PFNA),
            PFDA = min(PFDA), 
            PFUA = min(PFUA), .groups = 'keep') 

dsed.MAX<-d.sed  %>% 
  dplyr:::group_by(loc, Seasonality) %>% 
  summarize(
            PFHxS = max(PFHxS),
            PFOS = max(PFOS), 
            PFOA = max(PFOA), 
            PFNA = max(PFNA),
            PFDA = max(PFDA), 
            PFUA = max(PFUA), .groups = 'keep') 

# fish weight data, originally in g turn to kg
d.fish.w.wg$loc <- "WG"
d.fish.w.jba$Seasonality<-"na"
d.fish.w.jba$loc <- "JBA"
d.fish.w<-rbind(d.fish.w.jba, d.fish.w.wg)

# grams to kg
d.fish.w$mean_mass<-d.fish.w$mean_mass/1000
d.fish.w$min_mass<-d.fish.w$min_mass/1000
d.fish.w$max_mass<-d.fish.w$max_mass/1000

names(d.fish.w)[names(d.fish.w) == 'Sp'] <- 'SppAlias'


  
