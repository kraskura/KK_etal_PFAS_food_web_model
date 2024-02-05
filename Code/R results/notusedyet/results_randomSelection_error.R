# - spring vs summer 

library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse) 
library(readxl)
library(reshape)
library(rgl)
library(broom)
library(ggformat2)
library(cowplot)

# source all functions to run all model runs
source("./Code/PFAS_model_Sun_etal_classes_KK.R") # Sets up classes
source("./Code/PFAS_model_Sun_etal_steadyStateMod_KK.R") # steady state model, needed for 'runEcosystemModel' it is called in Bioaccumulation model 
source("./Code/PFAS_model_Sun_etal_BioaccumulationMod.R") # steady state model, needed for 'runEcosystemModel'
source("./Code/runEcosystemModels.R") 
options(error = traceback)

# Species, water, and sediment data
metadata<-as.data.frame(read_excel("./Data/Species_specific_info.xlsx", sheet = 1))
i<-c(1:9)
j<-c(10:ncol(metadata))
metadata[ , i] <- apply(metadata[ , i], 2,           
                    function(x) factor((x)))
metadata[ , j] <- apply(metadata[ , j], 2,            
                    function(x) as.numeric(x))

# Individual fish ID is now for sample ID, used to get individual means.  
metadata$SampleID<-gsub('-LIV','',metadata$SampleID) # liver
metadata$SampleID<-gsub('-GON','',metadata$SampleID) # gonad
metadata$SampleID<-gsub('-GN','',metadata$SampleID) # gonad
metadata$SampleID<-gsub('-MU','',metadata$SampleID) # muscle

metadata_means<-metadata %>% 
  dplyr:::group_by(SampleID, SppAlias) %>% 
  summarize(Weight_g = mean(Weight_g),
            PFHxS = mean(PFHxS),
            PFOS = mean(PFOS), 
            PFOA = mean(PFOA), 
            PFNA = mean(PFNA),
            PFDA = mean(PFDA), 
            PFUA = mean(PFUA)) 

# replace metadata with metadata means. 
for (i in 1:nrow(metadata)){
  r<-which(metadata_means$SampleID == metadata$SampleID[i])
  
  metadata$PFHxS[i] = metadata_means$PFHxS[r]
  metadata$PFOS[i] = metadata_means$PFOS[r]
  metadata$PFOA[i] = metadata_means$PFOA[r]
  metadata$PFNA[i] = metadata_means$PFNA[r]
  metadata$PFDA[i] = metadata_means$PFDA[r]
  metadata$PFUA[i] = metadata_means$PFUA[r]
}
  

# ********************************************************
WG_inputFiles_list = list(
    'numSpecies'= 4,
    'organismData'= read.csv('./Data/WG/organismData_BMFBCF_WG_fall.csv', row.names = "X"),
    'chemicalData'= read.csv('./Data/WG/chemicalData_full_WG_fall.csv', row.names = "chemicalParameter"),
    'oceanData'= read.csv('./Data/WG/oceanData_BCFBMF_WG_fall.csv', row.names = "oceanParameter"),
    'chemicalParams'= read.csv('./Data/WG/chemicalParameters_BCFBMF_WG.csv', row.names = "chemicalParameter"),
    'foodWebData'= read.csv('./Data/WG/foodWebTable_BCFBMF_WG_fish_fall.csv', row.names = "X"),
    'dietData' = read.csv('./Data/WG/dietData_median_WG_fall.csv', row.names = "Species"),  # changed to csv
    'min_dietData' = read.csv('./Data/WG/dietData_min_WG_fall.csv', row.names = "Species"),
    'max_dietData' = read.csv('./Data/WG/dietData_max_WG_fall.csv', row.names = "Species")
)


# spring
JBA_inputFiles_list = list(
    'numSpecies'= 9,
    'organismData'= read.csv('./Data/JBA/organismData_BMFBCF_JBA_spring.csv', row.names = "X"),
    'chemicalData'= read.csv('./Data/JBA/chemicalData_full_JBA_spring.csv', row.names = "chemicalParameter"),
    'oceanData'= read.csv('./Data/JBA/oceanData_BCFBMF_JBA_spring.csv', row.names = "oceanParameter"),
    'chemicalParams'= read.csv('./Data/JBA/chemicalParameters_BCFBMF_JBA.csv', row.names = "chemicalParameter"),
    'foodWebData'= read.csv('./Data/JBA/foodWebTable_BCFBMF_JBA_fish_springV2.csv', row.names = "X"),
    'dietData' = read.csv('./Data/JBA/dietData_median_JBA_spring.csv', row.names = "Species"),  # changed to csv
    'min_dietData' = read.csv('./Data/JBA/dietData_min_JBA_spring.csv', row.names = "Species"),
    'max_dietData' = read.csv('./Data/JBA/dietData_max_JBA_spring.csv', row.names = "Species")
)

# summer
JBA_inputFiles_list_summer = list(
    'numSpecies'= 9,
    'organismData'= read.csv('./Data/JBA/organismData_BMFBCF_JBA_summer.csv', row.names = "X"),
    'chemicalData'= read.csv('./Data/JBA/chemicalData_full_JBA_summer.csv', row.names = "chemicalParameter"),
    'oceanData'= read.csv('./Data/JBA/oceanData_BCFBMF_JBA_summer.csv', row.names = "oceanParameter"),
    'chemicalParams'= read.csv('./Data/JBA/chemicalParameters_BCFBMF_JBA.csv', row.names = "chemicalParameter"),
    'foodWebData'= read.csv('./Data/JBA/foodWebTable_BCFBMF_JBA_fish_summerV2.csv', row.names = "X"),
    'dietData' = read.csv('./Data/JBA/dietData_median_JBA_summer.csv', row.names = "Species"),  # changed to csv
    'min_dietData' = read.csv('./Data/JBA/dietData_min_JBA_summer.csv', row.names = "Species"),
    'max_dietData' = read.csv('./Data/JBA/dietData_max_JBA_summer.csv', row.names = "Species")
)


## 2. kRTable and settings - forced diets----
kRTable = read.csv('./Data/export_krTable.csv', row.names = "X", header = T)
colnames(kRTable)<-c("kr/kb", "PFAA", "chemID")

settings_obsDiet = new_Settings(
  chooseStudyType = 'default',
  chooseDiet = 'forced Munoz',
  chooseEd = 'Goeritz',
  chooseKoc = 'Koc_Munoz',
  chooseRenal ='off'
)
settings_modelDiet = new_Settings(
  chooseStudyType = 'default',
  chooseDiet = 'default',
  chooseEd = 'Goeritz',
  chooseKoc = 'Koc_Munoz',
  chooseRenal ='off'
)


# 3. Define changed parameters -----
# import data to make a random draw
temp_springJBA<-as.numeric(metadata[c(metadata$Site == "JBA" &
                             metadata$Season == "Spring" &
                             metadata$`Species/Sediment/Water` == "Water"), "Temperature"])
temp_sumJBA<-as.numeric(metadata[c(metadata$Site == "JBA" &
                             metadata$Season == "Summer" &
                             metadata$`Species/Sediment/Water` == "Water"), "Temperature"])
temp_WG<-as.numeric(metadata[c(metadata$Site == "WG" &
                             metadata$Season == "Fall" &
                             metadata$`Species/Sediment/Water` == "Water"), "Temperature"])

DO_springJBA<-as.numeric(metadata[c(metadata$Site == "JBA" &
                             metadata$Season == "Spring" &
                             metadata$`Species/Sediment/Water` == "Water"), "DO"])
DO_sumJBA<-as.numeric(metadata[c(metadata$Site == "JBA" &
                             metadata$Season == "Summer" &
                             metadata$`Species/Sediment/Water` == "Water"), "DO"])
DO_WG<-as.numeric(metadata[c(metadata$Site == "WG" &
                             metadata$Season == "Fall" &
                             metadata$`Species/Sediment/Water` == "Water"), "DO"])
# waters
Water_WG<-(metadata[c(metadata$Site == "WG" & metadata$`Species/Sediment/Water`== "Water"), c("Tissue/Location", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Water_JBAsp<-(metadata[c(metadata$Site == "JBA" & metadata$`Species/Sediment/Water`== "Water" & metadata$Season == "Spring"), c("Tissue/Location",'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Water_JBAsum<-(metadata[c(metadata$Site == "JBA" & metadata$`Species/Sediment/Water`== "Water" & metadata$Season == "Spring"), c("Tissue/Location", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])

# sediment
Sed_WG<-(metadata[c(metadata$Site == "WG" & metadata$`Species/Sediment/Water`== "Sediment"), c("Tissue/Location", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Sed_JBAsp<-(metadata[c(metadata$Site == "JBA" & metadata$`Species/Sediment/Water`== "Sediment" & metadata$Season == "Spring"), c("Tissue/Location",'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Sed_JBAsum<-(metadata[c(metadata$Site == "JBA" & metadata$`Species/Sediment/Water`== "Sediment" & metadata$Season == "Spring"), c("Tissue/Location", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])

# species subsets
# 1) WG
Pry_WG<-(metadata[c(metadata$Site == "WG" & metadata$SppAlias== "Pry"), c("SppAlias", "Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Bas_WG<-(metadata[c(metadata$Site == "WG" & metadata$SppAlias== "Bas"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Bgl_WG<-(metadata[c(metadata$Site == "WG" & metadata$SppAlias== "Bgl"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])

# 1) JBA  spring
Kil_JBAsp<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Kil" & metadata$Season == "Spring"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Pum_JBAsp<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Pum" & metadata$Season == "Spring"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Chu_JBAsp<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Chu" & metadata$Season == "Spring"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Dac_JBAsp<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Dac" & metadata$Season == "Spring"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Dar_JBAsp<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Dar" & metadata$Season == "Spring"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Min_JBAsp<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Min" & metadata$Season == "Spring"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Mad_JBAsp<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Mad" & metadata$Season == "Spring"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Swa_JBAsp<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Swa" & metadata$Season == "Spring"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
# 2) JBA  summer
Bas_JBAsum<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Bas" & metadata$Season == "Summer"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Pum_JBAsum<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Pum" & metadata$Season == "Summer"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Dac_JBAsum<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Dac" & metadata$Season == "Summer"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Dar_JBAsum<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Dar" & metadata$Season == "Summer"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Min_JBAsum<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Min" & metadata$Season == "Summer"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Mad_JBAsum<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Mad" & metadata$Season == "Summer"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Swa_JBAsum<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Swa" & metadata$Season == "Summer"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Fal_JBAsum<-(metadata[c(metadata$Site == "JBA" & metadata$SppAlias== "Fal" & metadata$Season == "Summer"), c("SppAlias","Tissue/Location", "Weight_g", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])



PFAA_cols <- c('PFHxS' = 'blueviolet', # 
             'PFOS' = 'cornflowerblue', # 'c'
             'PFOA' = 'orange',
             'PFNA' = 'violet',
             'PFDA' = 'deeppink',
             'PFUA' = 'limegreen')

# JBA: c(Kil, Chu, Dac, Dar, Min, Mad , Pum ,Swa)
# WG : c(Pry, Bgl, Bas)

remove(AllData)

# Loop for WG: one random dataset
species_list_WG<-c("Pry", "Bgl", "Bas")
species_list_JBAsp<-c("Kil", "Chu", "Dac", "Dar", "Min", "Mad", "Pum", "Swa")
species_list_JBAsum<-c("Dac", "Dar", "Min", "Fal", "Bas", "Mad", "Pum", "Swa")

species_data_list_WG<-list(Pry_WG, Bgl_WG, Bas_WG)
species_data_list_JBAsp<-list(Kil_JBAsp, Chu_JBAsp, Dac_JBAsp, Dar_JBAsp, Min_JBAsp, Mad_JBAsp, Pum_JBAsp, Swa_JBAsp)
species_data_list_JBAsum<-list(Dac_JBAsum, Dar_JBAsum, Min_JBAsum, Fal_JBAsum, Bas_JBAsum, Mad_JBAsum, Pum_JBAsum, Swa_JBAsum)

mytrophicLevel_WG = data.frame(SppAlias = c("Bas", "Bgl", "Phy", "Pry"), trophicLevel = c(3.69, 3.08, 1, 2.31))
mytrophicLevel_JBAsp = data.frame(SppAlias = c("Kil", "Chu", "Dac", "Dar", "Min", "Mad", "Pum", "Swa"),
                                  trophicLevel = c(3.4, 3.93, 2.62, 3.75, 2.92, 3.53, 3.83, 2.8)) # from Brown et al (Abbi's thesis)
mytrophicLevel_JBAsum = data.frame(SppAlias = c("Dac", "Dar", "Min", "Fal", "Bas", "Mad", "Pum", "Swa"),
                                   trophicLevel = c(2.62, 3.75, 2.92, 3.25, 3.8, 3.53, 3.83, 2.8))

# Swa :Trophic level (Ref. 69278):  2.8   ±0.28 se; based on food items. # https://www.fishbase.se/summary/2889
# Kil: Trophic level (Ref. 69278):  3.4   ±0.52 se; based on food items. # https://www.fishbase.se/summary/3191#:~:text=Trophic%20level%20(Ref.,se%3B%20based%20on%20food%20items.
# Bas: Trophic level (Ref. 69278):  3.8   ±0.4 se; based on diet studies.# https://www.fishbase.se/summary/micropterus-salmoides.html#:~:text=Trophic%20level%20(Ref.,(4.6%20%2D%206.1)%20years.


# function to get a random selection of parameters from observed data
rand_sel_inputFiles<-function(species,
                            inputFiles,
                            species_data,
                            temp_data,
                            sediment_data,
                            water_data){
  
  d0<-dplyr::sample_n(species_data, 1)
  t0<-sample(temp_data, 1)
  sed0<-dplyr::sample_n(sediment_data, 1)
  water0<-dplyr::sample_n(water_data, 1)
  inputFiles$organismData[species][1,1] <- d0[, "Weight_g"]/1000 # from kg to g
  inputFiles$oceanData["T", 1] <- t0
  
  inputFiles$dietData[species, ] <- d0[, c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")]
  inputFiles$chemicalData["C_s", 2:7]<-sed0[, c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")]
  inputFiles$chemicalData["C_WTO", 2:7]<-water0 [, c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")]/1000 # from ng/L to ng/ml
  
  return(list(inputFiles))

}
  
rand_ecosystem_runs<-function(n_iter,
                              ecosyst_ID,
                              mytrophicLevel,# Bas Bgl Phy Pry
                              species_list_ecosyst,
                              species_data_list_ecosyst, 
                              ecosyst_inputFiles_list,
                              temp_ecosyst,
                              Sed_ecosyst,
                              Water_ecosyst,
                              settings_Diet, # settings_obsDiet
                              ecosyst_min_dietData,
                              ecosyst_max_dietData, 
                              kRTabe = kRTable){
  
  for (i in 1:n_iter){ # loop for N random amount of datasets
    for(k in 1:length(species_list_ecosyst)){
      if(k == 1){
          inputFiles_ecosyst_rand<-rand_sel_inputFiles(species = species_list_ecosyst[k],
                                       species_data = species_data_list_ecosyst[[k]],
                                       inputFiles = ecosyst_inputFiles_list,
                                       temp_data = temp_ecosyst,
                                       sediment_data = Sed_ecosyst,
                                       water_data = Water_ecosyst)
      }else{
          inputFiles_ecosyst_rand<-rand_sel_inputFiles(species = species_list_ecosyst[k],
                                       species_data = species_data_list_ecosyst[[k]],
                                       inputFiles = inputFiles_ecosyst_rand[[1]],
                                       temp_data = temp_ecosyst,
                                       sediment_data = Sed_ecosyst,
                                       water_data = Water_ecosyst)
      }
    }
    inputFiles_ecosyst_rand<-inputFiles_ecosyst_rand[[1]]
    AllData.0<-runEcosystemModel(settings = settings_Diet,
                    inputFiles_list = inputFiles_ecosyst_rand,
                    # parameterList = c('C_B','C_WTO','C_s','k2','ke','kg','Total Elimination'),
                    parameterList = c("C_B","Phi", "mO" ,"C_WTO", "C_WDP" ,"Water", "C_s", "FeedRate",
                                      "Diet", "Gill_uptake", "Dietary_uptake","Gill_up%",
                                          "Diet_up%", "G_V", "G_D", "G_F", 
                                          "W_B", "Ew", "Ed", "k1",
                                          "k2", "kd", "ke", "kg",
                                          "kr_est", "Total Elimination", "kr_pct", "pKa", 
                                          "logDbw", "log Kow", "log Dmw", "log Dow", 
                                          "log Kpw", "D_BW", "D_MW", "D_OW",
                                          "K_PW", "pHi", "pHg", "K_GB",
                                          "nu_NB", "nu_LB", "nu_PB", "nu_OB", 
                                          "nu_WB", "epsilon_N", "epsilon_L", "epsilon_P",
                                          "epsilon_O", "Log_Koc", "Phi", "logC_B",
                                          "C_B_ngg", "logBAF", "BMF"),
                    PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                    dietData = inputFiles_ecosyst_rand$dietData,
                    kRTable = kRTable,
                    run_w_RenalElim = TRUE,
                    food_web_calc = TRUE,
                    min_diet = FALSE,
                    max_diet = FALSE,
                    min_dietData = ecosyst_min_dietData, 
                    max_dietData = ecosyst_max_dietData,
                    RunID = paste(ecosyst_ID, i, sep ="_"))
    if(i==1 & !exists("AllData")){
      AllData<-AllData.0
    }else{

      AllData<-rbind(AllData, AllData.0)
    }
    
    if(i == n_iter){
      AllData<-merge(AllData, mytrophicLevel, by = "SppAlias")
    }

  
  }
  
  AllData<-AllData[!c(AllData$PFAA=="PFUA" | AllData$PFAA=="PFDA"), ]
  AllData$C_B_re_ngg<- AllData$C_B_re / 1000 # was in ngkg 
  
  # stackbar dataset
  AllDataProp<-melt(AllData[, c("SppAlias", "PFAA", "Gill_up%", "Diet_up%", "runID")], id.vars=c("runID","SppAlias", "PFAA"))
  AllDataPropRate<-melt(AllData[, c("SppAlias", "PFAA", "Gill_uptake", "Dietary_uptake", "runID")], id.vars=c("runID","SppAlias", "PFAA"))
  AllDataPropRate2<-melt(AllData[, c("SppAlias", "PFAA", "k1", "k2", "runID")], id.vars=c("runID","SppAlias", "PFAA"))
  
  AllDataElProp<-melt(AllData[, c("SppAlias", "PFAA", "k2","ke", "kg","kr_est", "runID", "WB", "Total Elimination")], id.vars=c("runID","SppAlias", "PFAA", "Total Elimination", "WB"))
  
  AllData_sum_ValUp_mean<-AllData %>%
    dplyr::group_by(PFAA, SppAlias) %>% 
    summarize(Gill = mean(`Gill_uptake`),
              Diet = mean(`Dietary_uptake`)) %>% 
    pivot_longer(cols = c(Gill, Diet), names_to = "Uptake_route", values_to = "Uptake_g.kg.day") 
  AllData_sum_ValUp_sd<-AllData %>% 
    dplyr::group_by(PFAA, SppAlias) %>% 
    summarize(Gill = sd(`Gill_uptake`),
              Diet = sd(`Dietary_uptake`)) %>% 
    pivot_longer(cols = c(Gill, Diet), names_to = "Uptake_route", values_to = "Uptake_SD_g.kg.day") 
  AllData_meansVal <- merge(AllData_sum_ValUp_mean, AllData_sum_ValUp_sd, all.x = TRUE, by = c("PFAA", "SppAlias", "Uptake_route"))
  
  
  AllData_sum_PercUp_mean<-AllData %>% 
    dplyr::group_by(PFAA, SppAlias) %>% 
    summarize(Gill = mean(`Gill_up%`),
              Diet = mean(`Diet_up%`)) %>% 
    pivot_longer(cols = c(Gill, Diet), names_to = "Uptake_route", values_to = "Uptake_%g.kg.day") 
    AllData_sum_PercUp_mean
  AllData_sum_PercUp_sd<-AllData %>% 
    dplyr::group_by(PFAA, SppAlias) %>% 
    summarize(Gill = sd(`Gill_up%`),
              Diet = sd(`Diet_up%`)) %>% 
    pivot_longer(cols = c(Gill, Diet), names_to = "Uptake_route", values_to = "Uptake_SD_%g.kg.day") 

  AllData_means <- merge(AllData_sum_PercUp_mean, AllData_sum_PercUp_sd, all.x = TRUE, by = c("PFAA", "SppAlias", "Uptake_route"))
  
  # what factors contribute to higher difference
  AllData$diff_Obs_mod_logngkg_re<-AllData$Obs_logngkg - AllData$log_ngkg_re
  AllData$diff_Obs_mod_logngkg<-AllData$Obs_logngkg - AllData$log_ngkg
  
  # PERCENT 
  p1<-ggplot(AllData_means,
      aes(fill=Uptake_route, x=`Uptake_%g.kg.day`*100, y=SppAlias, color = Uptake_route)) + 
      geom_bar(stat="identity", position = "stack") +
      # geom_pointrange(aes(y = SppAlias,
      #                xmin=`Uptake_%g.kg.day`-`Uptake_SD_%g.kg.day`,
      #                xmax=`Uptake_%g.kg.day`+`Uptake_SD_%g.kg.day`), position = position_dodge(width = 1))+
      facet_wrap(.~PFAA, scales = "free")+
      scale_fill_manual(values = c("green3", "dodgerblue"))+
      scale_color_manual(values = c("green4", "darkblue"))+
      xlab("% Uptake")+
      ylab("Species")+
    theme_bw()+
    theme(legend.position = "top")
  
  
  p2<-ggplot(AllData[as.numeric(factor(AllData$runID)) == 1,],
      aes(x=Diet, y=SppAlias)) + 
      geom_bar(mapping = aes(x = C_B_re_ngg , y=SppAlias), 
               stat="identity", position = "stack", fill ="grey50", alpha = 0.5) +  # ngkg
      geom_bar(mapping = aes(x = C_B_re_ngg * `Diet_up%`, y=SppAlias), 
               stat="identity", position = "stack", fill ="grey50", alpha = 0.5) +  # ngkg
      # geom_bar(stat="identity", position = "stack", fill ="green4") + # g chemical/kg fish/day
       geom_bar(mapping = aes(x = Diet, y=SppAlias), # ng/g (or g/kg food) 
               stat="identity", position = "stack",fill ="green3", alpha = 0.5)+
     facet_wrap(.~PFAA, scales = "free")+
      ylab("Species")+
      xlab("PFAA levels [fish = ng/g, diet = ng/g]")+
    theme_bw()
  
  
  p3<-ggplot(AllData[as.numeric(factor(AllData$runID)) == 1,],
      aes(x=Gill_uptake, y=SppAlias)) + 
      geom_bar(mapping = aes(x = Water*1000, y=SppAlias), # ng/L (or g/kg food) 
               stat="identity", position = "stack",fill ="dodgerblue", alpha = 0.5)+
     geom_bar(mapping = aes(x = C_B_re_ngg , y=SppAlias), 
               stat="identity", position = "stack", fill ="grey50", alpha = 0.5) +  # ngkg
     geom_bar(mapping = aes(x = C_B_re_ngg * `Gill_up%`, y=SppAlias), 
               stat="identity", position = "stack", fill ="grey50", alpha = 0.5) +  # ngkg
      # geom_bar(stat="identity", position = "stack", fill ="green4") + # g chemical/kg fish/day
      facet_wrap(.~PFAA, scales = "free")+
      ylab("Species")+
      xlab("PFAA levels [fish = ng/g, water = ng/L]")+
    theme_bw()
  
  
  p4<-ggplot(AllData[as.numeric(factor(AllData$runID)) == 1,],
      aes(x=Gill_uptake, y=SppAlias)) + 
      geom_bar(mapping = aes(x = C_WTO*1000, y=SppAlias), # ng/L (or g/kg food) 
               stat="identity", position = position_nudge(y = 0.1), fill ="skyblue", width = 0.1, alpha = 1)+
     geom_bar(mapping = aes(x = C_B_re_ngg , y=SppAlias), 
               stat="identity",  position = position_nudge(y = 0), fill ="black", width = 0.1, alpha = 1) +  # ngkg
     geom_bar(mapping = aes(x = C_s, y=SppAlias), 
               stat="identity",   position = position_nudge(y = -0.1),  fill ="chocolate", width = 0.1, alpha = 1) +  # ngkg
     geom_bar(mapping = aes(x = Diet, y=SppAlias), 
               stat="identity",   position = position_nudge(y = -0.2),  fill ="green4", width = 0.1, alpha = 1) +  # ngkg
      # geom_bar(stat="identity", position = "stack", fill ="green4") + # g chemical/kg fish/day
      facet_wrap(.~PFAA, scales = "free")+
      ylab("Species")+
      xlab("PFAA: uptake (black), water (blue), diet (green), and sediment(brown)")+
      ggtitle(paste(ecosyst_ID))
    theme_bw()
  
  p5<-ggplot(data = AllData,
         aes(y = diff_Obs_mod_logngkg_re, x = as.numeric(trophicLevel),
             color = PFAA, alpha = temperature, group = runID))+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_point(pch=19)+
    geom_line(linewidth = 0.1)+
    facet_wrap(PFAA~., scale = "free")+
    scale_color_manual(values = PFAA_cols)+
    theme_classic()
  
  p6<-ggplot(data = AllData,
         aes(y = diff_Obs_mod_logngkg_re, x = `Diet_up%`,
             color = PFAA, alpha = temperature))+
    geom_hline(yintercept = 0, linetype = "dashed")+
    geom_point(pch=19)+
    geom_smooth(linewidth = 0.1, alpha=0.1)+
    facet_wrap(PFAA~SppAlias, scale = "free")+
    scale_color_manual(values = PFAA_cols)+
    theme_classic()
  
  p7<-ggplot(data = AllData,
         aes(y = logBAF, x = log10(BMF),
             color = PFAA, alpha = temperature))+
    geom_point(pch=19)+
    facet_wrap(.~SppAlias)+
    scale_color_manual(values = PFAA_cols)+
    theme_classic()
  
  p8<-ggplot(data = AllData,
         aes(y = logBAF/log10(BMF), x = Diet,
             color = PFAA, alpha = temperature))+
    geom_point(pch=19)+
    geom_smooth(linewidth = 0.1, alpha=0.1)+
    facet_wrap(PFAA~SppAlias, scales = "free")+
    scale_color_manual(values = PFAA_cols)+
    theme_classic()
  
  plot_save<-cowplot::plot_grid(p1, p4, p2, p3)

  ggsave(paste("./Figures/", ecosyst_ID, ".png",sep = ""), plot_save , width = 8, height = 8, units = "in")
  ggsave(paste("./Figures/", ecosyst_ID, "_trophicLev.png",sep = ""), p6 , width = 8, height = 8, units = "in")
  ggsave(paste("./Figures/", ecosyst_ID, "_logBAF_BMF.png",sep = ""), p7 , width = 5, height = 5, units = "in")

  return(list(AllData, p1, p2, p3, p4, p5, p6, p7, p8))
}


WG_observed<-rand_ecosystem_runs(n_iter=500,
                        ecosyst_ID = "WG",
                        mytrophicLevel = mytrophicLevel_WG,# Bas Bgl Phy Pry
                        species_list_ecosyst = species_list_WG,
                        species_data_list_ecosyst = species_data_list_WG, 
                        ecosyst_inputFiles_list = WG_inputFiles_list,
                        temp_ecosyst = temp_WG,
                        Sed_ecosyst = Sed_WG,
                        Water_ecosyst = Water_WG,
                        settings_Diet = settings_obsDiet, 
                        ecosyst_min_dietData = NULL,
                        ecosyst_max_dietData = NULL)

WG_modeled<-rand_ecosystem_runs(n_iter=500,
                        ecosyst_ID = "WG",
                        mytrophicLevel = mytrophicLevel_WG,# Bas Bgl Phy Pry
                        species_list_ecosyst = species_list_WG,
                        species_data_list_ecosyst = species_data_list_WG, 
                        ecosyst_inputFiles_list = WG_inputFiles_list,
                        temp_ecosyst = temp_WG,
                        Sed_ecosyst = Sed_WG,
                        Water_ecosyst = Water_WG,
                        settings_Diet = settings_modelDiet,
                        ecosyst_min_dietData = NULL,
                        ecosyst_max_dietData = NULL)

JBAsp_observed<-rand_ecosystem_runs(n_iter=50,
                        ecosyst_ID = "JBAsp",
                        mytrophicLevel = mytrophicLevel_JBAsp,
                        species_list_ecosyst = species_list_JBAsp,
                        species_data_list_ecosyst = species_data_list_JBAsp, 
                        ecosyst_inputFiles_list = JBA_inputFiles_list,
                        temp_ecosyst = temp_springJBA,
                        Sed_ecosyst = Sed_JBAsp,
                        Water_ecosyst = Water_JBAsp,
                        settings_Diet = settings_obsDiet, 
                        ecosyst_min_dietData = NULL,
                        ecosyst_max_dietData = NULL)

JBAsp_modeled<-rand_ecosystem_runs(n_iter=500,
                        ecosyst_ID = "JBAsp",
                        mytrophicLevel = mytrophicLevel_JBAsp,
                        species_list_ecosyst = species_list_JBAsp,
                        species_data_list_ecosyst = species_data_list_JBAsp, 
                        ecosyst_inputFiles_list = JBA_inputFiles_list,
                        temp_ecosyst = temp_springJBA,
                        Sed_ecosyst = Sed_JBAsp,
                        Water_ecosyst = Water_JBAsp,
                        settings_Diet = settings_modelDiet,
                        ecosyst_min_dietData = NULL,
                        ecosyst_max_dietData = NULL)

JBAsum_observed<-rand_ecosystem_runs(n_iter=500,
                        ecosyst_ID = "JBAsum",
                        mytrophicLevel = mytrophicLevel_JBAsum,
                        species_list_ecosyst = species_list_JBAsum,
                        species_data_list_ecosyst = species_data_list_JBAsum,
                        ecosyst_inputFiles_list = JBA_inputFiles_list_summer,
                        temp_ecosyst = temp_sumJBA,
                        Sed_ecosyst = Sed_JBAsum,
                        Water_ecosyst = Water_JBAsum,
                        settings_Diet = settings_obsDiet,
                        ecosyst_min_dietData = NULL,
                        ecosyst_max_dietData = NULL)

JBAsum_modeled<-rand_ecosystem_runs(n_iter=500,
                        ecosyst_ID = "JBAsum",
                        mytrophicLevel = mytrophicLevel_JBAsum,
                        species_list_ecosyst = species_list_JBAsum,
                        species_data_list_ecosyst = species_data_list_JBAsum,
                        ecosyst_inputFiles_list = JBA_inputFiles_list_summer,
                        temp_ecosyst = temp_sumJBA,
                        Sed_ecosyst = Sed_JBAsum,
                        Water_ecosyst = Water_JBAsum,
                        settings_Diet = settings_modelDiet,
                        ecosyst_min_dietData = NULL,
                        ecosyst_max_dietData = NULL)

JBAsum_observed[[6]]
JBAsum_modeled[[6]]

JBAsp_observed[[6]]
JBAsp_modeled[[6]]

WG_observed[[6]]
WG_modeled[[6]]


#### FIGURES for SERDP meeting (2023, july): -------
# 500 iterations

JBAsum_obsData<-JBAsum_observed[[1]]
JBAsum_modData<-JBAsum_modeled[[1]]
JBAsum_obsData<-JBAsum_obsData[!c(JBAsum_obsData$Obs_ngg == 0), ]
JBAsum_modData<-JBAsum_modData[!c(JBAsum_modData$Obs_ngg == 0), ]
write.csv(file = "./Data/poster_JBAsum_obsData_jun25.csv", x = JBAsum_obsData)
write.csv(file = "./Data/poster_JBAsum_modData_jun25.csv", x = JBAsum_modData)

JBAsp_obsData<-JBAsp_observed[[1]]
JBAsp_modData<-JBAsp_modeled[[1]]
JBAsp_obsData<-JBAsp_obsData[!c(JBAsp_obsData$Obs_ngg == 0), ]
JBAsp_modData<-JBAsp_modData[!c(JBAsp_modData$Obs_ngg == 0), ]
write.csv(file = "./Data/poster_JBAsp_obsData_jun25.csv", x = JBAsp_obsData)
write.csv(file = "./Data/poster_JBAsp_modData_jun25.csv", x = JBAsp_modData)

WG_obsData<-WG_observed[[1]]
WG_modData<-WG_modeled[[1]]
WG_obsData<-WG_obsData[!c(WG_obsData$Obs_ngg == 0), ]
WG_modData<-WG_modData[!c(WG_modData$Obs_ngg == 0), ]
write.csv(file = "./Data/poster_WG_obsData_jun25.csv", x = WG_obsData)
write.csv(file = "./Data/poster_WG_modData_jun25.csv", x = WG_modData)


# stackbar dataset
JBAsum_obsData_m<-melt(JBAsum_obsData[, c("SppAlias", "PFAA", "Gill_up%", "Diet_up%", "runID")], id.vars=c("runID","SppAlias", "PFAA"))
AllData_sum_ValUp_mean<-JBAsum_obsData %>%
dplyr::group_by(PFAA, SppAlias) %>%
summarize(Gill = mean(`Gill_uptake`),
          Diet = mean(`Dietary_uptake`)) %>%
pivot_longer(cols = c(Gill, Diet), names_to = "Uptake_route", values_to = "Uptake_g.kg.day")
AllData_sum_ValUp_sd<-JBAsum_obsData %>%
dplyr::group_by(PFAA, SppAlias) %>%
summarize(Gill = sd(`Gill_uptake`),
          Diet = sd(`Dietary_uptake`)) %>%
pivot_longer(cols = c(Gill, Diet), names_to = "Uptake_route", values_to = "Uptake_SD_g.kg.day")
AllData_meansVal <- merge(AllData_sum_ValUp_mean, AllData_sum_ValUp_sd, all.x = TRUE, by = c("PFAA", "SppAlias", "Uptake_route"))

# 
# AllData_sum_PercUp_mean<-AllData %>% 
# dplyr::group_by(PFAA, SppAlias) %>% 
# summarize(Gill = mean(`Gill_up%`),
#           Diet = mean(`Diet_up%`)) %>% 
# pivot_longer(cols = c(Gill, Diet), names_to = "Uptake_route", values_to = "Uptake_%g.kg.day") 
# AllData_sum_PercUp_mean
# AllData_sum_PercUp_sd<-AllData %>% 
# dplyr::group_by(PFAA, SppAlias) %>% 
# summarize(Gill = sd(`Gill_up%`),
#           Diet = sd(`Diet_up%`)) %>% 
# pivot_longer(cols = c(Gill, Diet), names_to = "Uptake_route", values_to = "Uptake_SD_%g.kg.day") 
# 
# AllData_means <- merge(AllData_sum_PercUp_mean, AllData_sum_PercUp_sd, all.x = TRUE, by = c("PFAA", "SppAlias", "Uptake_route"))
#   # what factors contribute to higher difference
#   AllData$diff_Obs_mod_logngkg_re<-AllData$Obs_logngkg - AllData$log_ngkg_re
#   AllData$diff_Obs_mod_logngkg<-AllData$Obs_logngkg - AllData$log_ngkg

p_lines <- data.frame('x' = c(1, 6),'y' = c(1, 6))
p_lines$y1 <- p_lines$x-1
p_lines$y2 <- p_lines$x+1
p_lines$y3 <- p_lines$x-log10(2)
p_lines$y4 <- p_lines$x+log10(2)

  # PERCENT 
p1<-ggplot(JBAsum_obsData_m[JBAsum_obsData_m$runID=="JBAsum_1", ],
      aes(fill=variable, x=value*100, y=SppAlias, color = variable)) + 
    geom_bar(stat="identity", position = "stack") +
    # geom_pointrange(aes(y = SppAlias,
    #                xmin=`Uptake_%g.kg.day`-`Uptake_SD_%g.kg.day`,
    #                xmax=`Uptake_%g.kg.day`+`Uptake_SD_%g.kg.day`), position = position_dodge(width = 1))+
    facet_wrap(.~PFAA, scales = "free")+
    scale_fill_manual(values = c("green3", "dodgerblue"),  labels=c("Gill uptake", "Diet uptake"))+
    scale_color_manual(values = c("green4", "darkblue"), labels=c("Gill uptake", "Diet uptake"))+
  theme_bw()
p1<-ggformat(p1, y_title = "Species", x_title = "% Uptake", print = F)
p1<-p1 + theme(legend.position = "top", legend.title = element_blank(), legend.text=element_text(size=15))
ggsave(filename = "./Figures/poster_Fig1_JBAsummer_observedDiet.png", plot = p1, width = 6, height = 7)


##### JBA SUMMER ------
JBAsum_modData_p<-
    ggplot(data = JBAsum_modData,
         aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, group = runID))+
    theme_classic()+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 2, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 2, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 1, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 1, linewidth = 0.4, color = "grey30")+
    ggalt::geom_encircle(mapping = aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, fill = PFAA, group = PFAA),
                         alpha=0.3, s_shape=1, expand=0, size=1)+
    scale_color_manual(values = PFAA_cols)+
    scale_fill_manual(values = PFAA_cols)+
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black", size = 1.5)+
    # theme(legend.position = "none")+
    geom_point(pch=19, size=1, alpha = 0.7)+
    ylim(0, 7)+
    xlim(0, 7)
JBAsum_modData_p<-ggformat(JBAsum_modData_p, x_title = 'Oberved log(ng/kg)', y_title = 'Modeled log(ng/kg)', size_text = 23, print = T)
JBAsum_modData_p<-JBAsum_modData_p+theme(legend.position = "none")
ggsave(filename = "./Figures/poster_Fig2a_JBAsummer_modeledDiet.png", plot = JBAsum_modData_p, width = 6, height = 6, dpi = 300)

JBAsum_obsData_p<-
    ggplot(data = JBAsum_obsData,
         aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, group = runID))+
    theme_classic()+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 2, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 2, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 1, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 1, linewidth = 0.4, color = "grey30")+
    ggalt::geom_encircle(mapping = aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, fill = PFAA, group = PFAA),
                         alpha=0.3, s_shape=1, expand=0, size=1)+
    scale_color_manual(values = PFAA_cols)+
    scale_fill_manual(values = PFAA_cols)+
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black", size = 1.5)+
    # theme(legend.position = "none")+
    geom_point(pch=19, size=1, alpha = 0.7)+
    ylim(0, 7)+
    xlim(0, 7)
JBAsum_obsData_p<-ggformat(JBAsum_obsData_p, x_title = 'Oberved log(ng/kg)', y_title = 'Modeled log(ng/kg)', size_text = 23, print = T)
JBAsum_obsData_p<-JBAsum_obsData_p+theme(legend.position = "none")
ggsave(filename = "./Figures/poster_Fig2b_JBAsummer_observedDiet.png", plot = JBAsum_obsData_p, width = 6, height = 6, dpi = 300)

##### JBA spring -------
JBAsp_modData_p<-
    ggplot(data = JBAsp_modData,
         aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, group = runID))+
    theme_classic()+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 2, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 2, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 1, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 1, linewidth = 0.4, color = "grey30")+
    ggalt::geom_encircle(mapping = aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, fill = PFAA, group = PFAA),
                         alpha=0.3, s_shape=1, expand=0, size=1)+
    scale_color_manual(values = PFAA_cols)+
    scale_fill_manual(values = PFAA_cols)+
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black", size = 1.5)+
    # theme(legend.position = "none")+
    geom_point(pch=19, size=1, alpha = 0.7)+
    ylim(0, 7)+
    xlim(0, 7)
JBAsp_modData_p<-ggformat(JBAsp_modData_p, x_title = 'Oberved log(ng/kg)', y_title = 'Modeled log(ng/kg)', size_text = 23, print = T)
JBAsp_modData_p<-JBAsp_modData_p+theme(legend.position = "none")
ggsave(filename = "./Figures/poster_Fig3a_JBAspring_modeledDiet.png", plot = JBAsp_modData_p, width = 6, height = 6, dpi = 300)

JBAsp_obsData_p<-
    ggplot(data = JBAsp_obsData,
         aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, group = runID))+
    theme_classic()+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 2, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 2, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 1, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 1, linewidth = 0.4, color = "grey30")+
    ggalt::geom_encircle(mapping = aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, fill = PFAA, group = PFAA),
                         alpha=0.3, s_shape=1, expand=0, size=1)+
    scale_color_manual(values = PFAA_cols)+
    scale_fill_manual(values = PFAA_cols)+
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black", size = 1.5)+
    # theme(legend.position = "none")+
    geom_point(pch=19, size=1, alpha = 0.7)+
    ylim(0, 7)+
    xlim(0, 7)
JBAsp_obsData_p<-ggformat(JBAsp_obsData_p, x_title = 'Oberved log(ng/kg)', y_title = 'Modeled log(ng/kg)', size_text = 23, print = T)
JBAsp_obsData_p<-JBAsp_obsData_p+theme(legend.position = "none")
ggsave(filename = "./Figures/poster_Fig3b_JBAspring_observedDiet.png", plot = JBAsp_obsData_p, width = 6, height = 6, dpi = 300)

##### JBA spring -------
WG_modData_p<-
    ggplot(data = WG_modData,
         aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, group = runID))+
    theme_classic()+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 2, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 2, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 1, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 1, linewidth = 0.4, color = "grey30")+
    ggalt::geom_encircle(mapping = aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, fill = PFAA, group = PFAA),
                         alpha=0.3, s_shape=1, expand=0, size=1)+
    scale_color_manual(values = PFAA_cols)+
    scale_fill_manual(values = PFAA_cols)+
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black", size = 1.5)+
    # theme(legend.position = "none")+
    geom_point(pch=19, size=1, alpha = 0.7)+
    ylim(0, 7)+
    xlim(0, 7)
WG_modData_p<-ggformat(WG_modData_p, x_title = 'Oberved log(ng/kg)', y_title = 'Modeled log(ng/kg)', size_text = 23, print = T)
WG_modData_p<-WG_modData_p+theme(legend.position = "none")
ggsave(filename = "./Figures/poster_Fig4a_WGfall_modeledDiet.png", plot = WG_modData_p, width = 6, height = 6, dpi = 300)

WG_obsData_p<-
    ggplot(data = WG_obsData,
         aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, group = runID))+
    theme_classic()+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 2, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 2, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 1, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 1, linewidth = 0.4, color = "grey30")+
    ggalt::geom_encircle(mapping = aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, fill = PFAA, group = PFAA),
                         alpha=0.3, s_shape=1, expand=0, size=1)+
    scale_color_manual(values = PFAA_cols)+
    scale_fill_manual(values = PFAA_cols)+
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black", size = 1.5)+
    # theme(legend.position = "none")+
    geom_point(pch=19, size=1, alpha = 0.7)+
    ylim(0, 7)+
    xlim(0, 7)
WG_obsData_p<-ggformat(WG_obsData_p, x_title = 'Oberved log(ng/kg)', y_title = 'Modeled log(ng/kg)', size_text = 23, print = T)
WG_obsData_p<-WG_obsData_p+theme(legend.position = "none")
ggsave(filename = "./Figures/poster_Fig4b_WGfall_observedDiet.png", plot = WG_obsData_p, width = 6, height = 6, dpi = 300)


##### legend only ------

legend_plot<-ggplot(data = WG_obsData,
         aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, group = runID))+
    theme_classic()+
    geom_point()+
    geom_line()+
    ggalt::geom_encircle(mapping = aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, fill = PFAA, group = PFAA),
                         alpha=0.3, s_shape=1, expand=0, size=1)+
    scale_color_manual(values = PFAA_cols)+
    scale_fill_manual(values = PFAA_cols)+
    theme(legend.title = element_blank(), legend.text=element_text(size=15))


legend <- cowplot::get_legend(legend_plot)
grid.newpage()
grid.draw(legend)

ggsave(filename = "./Figures/poster_Fig2-3-4_legend.png", plot = legend , width = 3, height = 3, dpi = 300)


ggplot(WG_obsData[c(WG_obsData$PFAA=="PFOS" & WG_obsData$SppAlias =="Bas"),],
               aes(y=`Dietary_uptake`/`Gill_uptake`, x = C_WTO , fill = C_s))+
  geom_point(pch=21,size= 3,  color ="black")+
  scale_fill_viridis_c(option = "C")+
  theme_bw()+
  xlab("PFOS in sediment (ng/g)")+
  ylab("PFOS ratio Diet:Water uptake")+
  ggtitle("WG bass only (field Sed, 10 ng/g water)")

ggplot(JBAsp_obsData,
               aes(alpha =`Diet_up%`*100, size = Dietary_uptake,
                   fill = PFAA,
                   y = C_WTO,
                   x = C_s))+
  geom_point( shape = 21,  color ="black")+
  scale_fill_manual(values = PFAA_cols)+ 
  guides(size = guide_legend("Dietary uptake (%)"), 
         alpha = guide_legend("Dietary uptake (ng/g)"))+
  ggalt::geom_encircle(data = JBAsp_obsData[JBAsp_obsData$PFAA == "PFOS",],
                   mapping =aes(alpha =`Diet_up%`*100, size = Dietary_uptake,
                   fill = PFAA,
                   y = C_WTO,
                   x = C_s, group = PFAA), alpha=0.5, s_shape=1, expand=0, size=1)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.5)+
  theme_bw()+
  ylim(0, 2.3)+
  xlim(0, 2.3)+
  facet_wrap(.~SppAlias)+
  xlab("PFAA in sediment (ng/g)")+
  ylab("PFAA in water (ng/mL)")+
  ggtitle("JBA spring (500 iterations)")
ggsave(filename = "./Figures/poster_JBAspring_sed_water.png", width = 8, height = 7)


library(car)
scatter3d(y=JBAsum_obsData[c(JBAsum_obsData$PFAA=="PFOS" & JBAsum_obsData$SppAlias =="Bas"),"Diet_up%"]*100,
      x=JBAsum_obsData[c(JBAsum_obsData$PFAA=="PFOS" & JBAsum_obsData$SppAlias =="Bas"),"C_s"], 
      z=JBAsum_obsData[c(JBAsum_obsData$PFAA=="PFOS" & JBAsum_obsData$SppAlias =="Bas"),"C_WTO"],
      point.col = "blue", 
        surface=T,
        ylab = "% dietary uptake",
        xlab = "PFOS in sediment (ng/g)",
        zlab = "PFOS in water (ng/mL)")

# ****************************************
# ****************************************
#### TEST HIGH SEDIMENT LOW WATER  --------------
# water - natural 
# make sure this is in ng/mL
Water_WG<-(metadata[c(metadata$Site == "WG" & metadata$`Species/Sediment/Water`== "Water"), c("Tissue/Location", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])/1000
Water_JBAsp<-(metadata[c(metadata$Site == "JBA" & metadata$`Species/Sediment/Water`== "Water" & metadata$Season == "Spring"), c("Tissue/Location",'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])/1000
Water_JBAsum<-(metadata[c(metadata$Site == "JBA" & metadata$`Species/Sediment/Water`== "Water" & metadata$Season == "Spring"), c("Tissue/Location", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])/1000

# water - artificial 
# sediment - artificial 
Water_WG$PFOS<-c(10, 10, 10, 10)
Water_JBAsp$PFOS<-c(rep(10, 24))
Water_JBAsum$PFOS<-c(rep(10, 24))

# water - artificial2
# sediment - natural 
Sed_WG<-(metadata[c(metadata$Site == "WG" & metadata$`Species/Sediment/Water`== "Sediment"), c("Tissue/Location", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Sed_JBAsp<-(metadata[c(metadata$Site == "JBA" & metadata$`Species/Sediment/Water`== "Sediment" & metadata$Season == "Spring"), c("Tissue/Location",'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])
Sed_JBAsum<-(metadata[c(metadata$Site == "JBA" & metadata$`Species/Sediment/Water`== "Sediment" & metadata$Season == "Spring"), c("Tissue/Location", 'PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')])

# # sediment - artificial 
# Sed_WG$PFOS<-rep(max(Sed_WG$PFOS), 4)
# Sed_JBAsp$PFOS<-rep(max(Sed_JBAsp$PFOS), 7)
# Sed_JBAsum$PFOS<-rep(max(Sed_JBAsum$PFOS), 7)

# artificial 2
# make sure this is in ng/mL
# Sed_WG$PFOS<-c(1, 100, 300, 400) / 1000
# Sed_JBAsp$PFOS<-c(1, 10, 200, 400, 600, 800, 900 ) / 1000
# Sed_JBAsum$PFOS<-c(1, 10, 200, 400, 600, 800, 900 ) / 1000

# Focus on JBA summer
JBAsum_modeled_SED<-rand_ecosystem_runs(n_iter=20,
                        ecosyst_ID = "JBAsum_SED",
                        mytrophicLevel = mytrophicLevel_JBAsum,
                        species_list_ecosyst = species_list_JBAsum,
                        species_data_list_ecosyst = species_data_list_JBAsum,
                        ecosyst_inputFiles_list = JBA_inputFiles_list_summer,
                        temp_ecosyst = temp_sumJBA,
                        Sed_ecosyst = Sed_JBAsum,
                        Water_ecosyst = Water_JBAsum,
                        settings_Diet = settings_modelDiet,
                        ecosyst_min_dietData = NULL,
                        ecosyst_max_dietData = NULL)

JBAsp_modeled_SED<-rand_ecosystem_runs(n_iter=20,
                        ecosyst_ID = "JBAsp_SED",
                        mytrophicLevel = mytrophicLevel_JBAsp,
                        species_list_ecosyst = species_list_JBAsp,
                        species_data_list_ecosyst = species_data_list_JBAsp,
                        ecosyst_inputFiles_list = JBA_inputFiles_list,
                        temp_ecosyst = temp_springJBA,
                        Sed_ecosyst = Sed_JBAsp,
                        Water_ecosyst = Water_JBAsp,
                        settings_Diet = settings_modelDiet,
                        ecosyst_min_dietData = NULL,
                        ecosyst_max_dietData = NULL)

WG_modeled_SED<-rand_ecosystem_runs(n_iter=20,
                        ecosyst_ID = "WG_SED",
                        mytrophicLevel = mytrophicLevel_WG,# Bas Bgl Phy Pry
                        species_list_ecosyst = species_list_WG,
                        species_data_list_ecosyst = species_data_list_WG,
                        ecosyst_inputFiles_list = WG_inputFiles_list,
                        temp_ecosyst = temp_WG,
                        Sed_ecosyst = Sed_WG,
                        Water_ecosyst = Water_WG,
                        settings_Diet = settings_modelDiet,
                        ecosyst_min_dietData = NULL,
                        ecosyst_max_dietData = NULL)

#  ****************************************
#  ****************************************

# generate a figure for these species (or just bass) in which you drop water concentrations way down
# but keep sediment constant at the “field value”? 

data_JBAsum<-JBAsum_modeled_SED[[1]]
p_sed1<-ggplot(data_JBAsum[data_JBAsum$PFAA=="PFOS",], aes(y=`Diet_up%`*100, x = C_s, color = SppAlias))+
  geom_point()+
  facet_grid(.~PFAA, scales = "free")+
  scale_color_viridis_d()+
  geom_smooth(method = "lm", se = FALSE)+
  # stat_smooth(method="lm", se=TRUE, fill=NA,
  #               formula=y ~ poly(x, 4, raw=TRUE))+
  ylim(0,100)+
  theme_bw()+
  xlab("PFOS in sediment (ng/g)")+
  ylab("PFOS dietary uptake (%)")+
  ggtitle("JBA summer (field Sed, 10 ng/g water)")
ggsave(filename = "./Figures/poster_sediment_PFOS_bioaccum_JBAsummer.png", width = 5, height = 5)
  

data_JBAsp<-JBAsp_modeled_SED[[1]]
p_sed2<-ggplot(data_JBAsp[data_JBAsp$PFAA=="PFOS",], aes(y=`Diet_up%`*100, x = C_s, color = SppAlias))+
  geom_point()+
  facet_grid(.~PFAA, scales = "free")+
  scale_color_viridis_d()+
  geom_smooth(method = "glm", se = FALSE)+
  ylim(0,100)+
  theme_bw()+
  xlab("PFOS in sediment (ng/g)")+
  ylab("PFOS dietary uptake (%)")+
  ggtitle("JBA spring (field Sed, 10 ng/g water)")
ggsave(filename = "./Figures/poster_sediment_PFOS_bioaccum_JBAspring.png", width = 5, height = 5)
 
data_WG<-WG_modeled_SED[[1]]
p_sed3<-ggplot(data_WG[c(data_WG$PFAA=="PFOS" & !c(data_WG$SppAlias =="Phy")),], aes(y=`Diet_up%`*100, x = C_s, color = SppAlias))+
  geom_point()+
  facet_grid(.~PFAA, scales = "free")+
  scale_color_viridis_d()+
  geom_smooth(method = "loess", se = FALSE)+
  ylim(0,100)+
  theme_bw()+
  xlab("PFOS in sediment (ng/g)")+
  ylab("PFOS dietary uptake (%)")+
  ggtitle("WG spring (field Sed, 10 ng/g water)")
ggsave(filename = "./Figures/poster_sediment_PFOS_bioaccum_WG.png", width = 5, height = 5)
 
# bass only 

data_JBAsum<-JBAsum_modeled_SED[[1]]
p_sedb<-ggplot(data_JBAsum[c(data_JBAsum$PFAA=="PFOS" & data_JBAsum$SppAlias =="Bas"),], aes(y=`Diet_up%`*100, x = C_s, color = SppAlias))+
  geom_point()+
  facet_grid(.~PFAA, scales = "free")+
  scale_color_viridis_d()+
  ylim(0,100)+
  geom_smooth(method = "lm", se = FALSE)+
  # stat_smooth(method="lm", se=TRUE, fill=NA,
  #               formula=y ~ poly(x, 4, raw=TRUE))+
  theme_bw()+
  xlab("PFOS in sediment (ng/g)")+
  ylab("PFOS dietary uptake (%)")+
  ggtitle("JBA summer (field Sed, 10 ng/g water)")
ggsave(filename = "./Figures/poster_sediment_PFOS_bioaccum_JBAsummer_Bass.png", width = 5, height = 5)
  

data_WG<-WG_modeled_SED[[1]]
p_sed3<-ggplot(data_WG[c(data_WG$PFAA=="PFOS" & data_WG$SppAlias =="Bas"),], aes(y=`Diet_up%`*100, x = C_s, color = SppAlias))+
  geom_point()+
  facet_grid(.~PFAA, scales = "free")+
  scale_color_viridis_d()+
  geom_smooth(method = "loess", se = FALSE)+
  ylim(0,100)+
  theme_bw()+
  xlab("PFOS in sediment (ng/g)")+
  ylab("PFOS dietary uptake (%)")+
  ggtitle("WG spring (field Sed, 10 ng/g water)")
ggsave(filename = "./Figures/poster_sediment_PFOS_bioaccum_WG_bass.png", width = 5, height = 5)
 

# generate a figure for just bass in which you somehow plot relative contribution of diet and waterborne 
# over a range of both of sediment in water 
ggplot(data_WG[c(data_WG$PFAA=="PFOS" & data_WG$SppAlias =="Bas"),],
               aes(y=`Dietary_uptake`/`Gill_uptake`, x = C_s , fill = C_s))+
  geom_point(pch=21,size= 3,  color ="black")+
  scale_fill_viridis_c(option = "C")+
  theme_bw()+
  xlab("PFOS in sediment (ng/g)")+
  ylab("PFOS ratio Diet:Water uptake")+
  ggtitle("WG bass only (field Sed, 10 ng/g water)")
# ggsave(filename = "./Figures/poster_sediment_test.png", width = 5, height = 5)





# 4. post run stats and error estimates --------
# 1. Quantify variation in the between the estimated values and observed values 

rand_ecosystem_runs_analysis<-function(data, obs_metadata, plot_title, Site, Season){
  
  metadata_l<-melt(obs_metadata,
      id.vars=c("Weight_g", "SppAlias", "Site", "Season"),
      measure.vars=c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA"),
  )
  colnames(metadata_l)<-c("Weight_g", "SppAlias", "Site", "Season", "PFAA", "value" )
  
  obs<-metadata_l %>% 
    filter(Site == Site, Season == Season) %>% 
    group_by(PFAA, SppAlias) %>% 
    summarise(mean = mean(value), 
              min = min(value),
              max = max(value),
              sd = sd(value),
              CV = (sd/mean), 
              var = var(value),
              n = n()) 
  
  mod<-data %>% 
    group_by(PFAA, SppAlias) %>% 
    summarise(mean = mean(C_B_re_ngg), 
              min = min(C_B_re_ngg),
              max = max(C_B_re_ngg),
              sd = sd(C_B_re_ngg),
              CV = (sd/mean), 
              var = var(C_B_re),
              n = n())
  
  var_stats<-merge(obs, mod, by = c("SppAlias", "PFAA"), all.y = T)
  
  # linear regression estimates  --------
  data[sapply(data, is.infinite)] <- NA
  
  modeled_reg_vals<-data %>%
    nest_by(runID) %>%
    mutate(mod = list(lm(log_ngkg_re ~ Obs_logngkg, data = data))) %>% 
    summarise(augment(mod)) %>% 
    as.data.frame()
  modeled_reg_vals<-merge(modeled_reg_vals, data[, c("PFAA","SppAlias", "runID", "Obs_logngkg")])
  modeled_reg_vals$diff_Obs_mod_logngkg<-modeled_reg_vals$log_ngkg_re - modeled_reg_vals$Obs_logngkg
  
  modeled_reg_comp<-data %>%
    nest_by(runID) %>%
    mutate(mod = list(lm(log_ngkg_re ~ Obs_logngkg, data = data))) %>% 
    summarise(glance(mod)) %>% 
    as.data.frame()
  
  modeled_reg_coefs<-data %>%
    nest_by(runID) %>%
    mutate(mod = list(lm(log_ngkg_re ~ Obs_logngkg, data = data))) %>% 
    summarise(tidy(mod)) %>% 
    as.data.frame()
  
  # Stats on the estimated parameters
  sum_est <- modeled_reg_coefs %>% 
    group_by(term) %>% 
    summarise(mean_val = mean(estimate), 
              sd_val = sd(estimate))
  
  ### plots -----
  p_cv<-ggplot(var_stats, aes(CV.x * 100, CV.y * 100, color = PFAA, shape = SppAlias))+
    geom_point()+
    scale_color_manual(values = PFAA_cols, guide = "none")+
    ylim(0, 300)+
    xlim(0, 300)+
    ylab("CV % (predicted fit data)")+
    xlab("CV % (observed data)")+
    theme_bw()+
    theme(legend.position = "top")
  
  p_mean<-ggplot(var_stats, aes(mean.x, mean.y, color = PFAA, shape = SppAlias))+
    geom_point()+
    geom_pointrange(mapping = aes(x=mean.x, y=mean.y, ymin=min.y, ymax=max.y))+
    geom_pointrange(mapping = aes(x=mean.x, y=mean.y, xmin=min.x, xmax=max.x))+
    scale_color_manual(values = PFAA_cols, guide = "none")+
    facet_wrap(PFAA~., scales = "free")+
    # ylim(0, 300)+
    # xlim(0, 300)+
    ylab("mean +/- range (predicted fit data)")+
    xlab("mean +/- range (observed data)")+
    theme_bw()+
    theme(legend.position = "none")
  
  # residuals 
  p_resid<-ggplot(modeled_reg_vals, aes(x = diff_Obs_mod_logngkg, fill = PFAA))+
    geom_histogram()+
    geom_vline(xintercept = 0, color = "red")+
    scale_color_manual(values = PFAA_cols)+
    scale_fill_manual(values = PFAA_cols)+
    theme_bw()+
    facet_grid(.~PFAA)+
    ylab("Count")+
    xlab("Residual difference (run regr.)")+
    xlim(-1.5, 1.5)+
    theme(legend.position = "none")
  
  p_resid1<-ggplot(modeled_reg_vals, aes(x = .resid))+
    geom_histogram(color = "grey50", fill = "grey50")+
    geom_vline(xintercept = 0, color = "red")+
    theme_bw()+
    facet_grid(.~PFAA)+
    ylab("Count")+
    xlab("Residual difference (LOI)")+
    xlim(-1.5, 1.5)+
    theme(legend.position = "none")
  
  # slopes and intercepts
  pred_int<-ggplot(modeled_reg_coefs[(modeled_reg_coefs$term == "(Intercept)"), ], aes(x = estimate, group = term))+
    geom_histogram(alpha = 0.5, color = "grey")+
    geom_vline(xintercept = sum_est$mean_val[1], size = 1, color = "black")+
    geom_vline(xintercept = sum_est$mean_val[1]+ sum_est$sd_val[1],
                size = 0.3, linetype = "dashed", color = "black")+
    geom_vline(xintercept = sum_est$mean_val[1]- sum_est$sd_val[1],
                size = 0.3, linetype = "dashed", color = "black")+
    geom_vline(xintercept = 0, color = "red")+
    theme_bw()+
    ylab("Count")+
    xlab("(Intercept)")+
    xlim(-0.1, 3)
  
  # slopes and intercepts
  pred_b<-ggplot(modeled_reg_coefs[(modeled_reg_coefs$term == "Obs_logngkg"), ], aes(x = estimate, group = term))+
    geom_histogram(alpha = 0.5, color = "grey")+
    geom_vline(xintercept = sum_est$mean_val[2], size = 1, color = "black")+
    geom_vline(xintercept = sum_est$mean_val[2]+ sum_est$sd_val[2],
                size = 0.3, linetype = "dashed", color = "black")+
    geom_vline(xintercept = sum_est$mean_val[2]- sum_est$sd_val[2],
                size = 0.3, linetype = "dashed", color = "black")+
    geom_vline(xintercept = 1, color = "red")+
    theme_bw()+
    ylab("Count")+
    xlab("Slope [Obs_logngkg]")+
    xlim(0, 1.1)
  
    pred_b<-ggplot(modeled_reg_coefs[(modeled_reg_coefs$term == "Obs_logngkg"), ], aes(x = estimate, group = term))+
      geom_histogram(alpha = 0.5, color = "grey")+
      geom_vline(xintercept = sum_est$mean_val[2], size = 1, color = "black")+
      geom_vline(xintercept = sum_est$mean_val[2]+ sum_est$sd_val[2],
                  size = 0.3, linetype = "dashed", color = "black")+
      geom_vline(xintercept = sum_est$mean_val[2]- sum_est$sd_val[2],
                  size = 0.3, linetype = "dashed", color = "black")+
      geom_vline(xintercept = 1, color = "red")+
      theme_bw()+
      ylab("Count")+
      xlab("Slope [Obs_logngkg]")+
      xlim(0, 1.1)

  p_lines <- data.frame('x' = c(1, 6),'y' = c(1, 6))
  p_lines$y1 <- p_lines$x-1
  p_lines$y2 <- p_lines$x+1
  p_lines$y3 <- p_lines$x-log10(2)
  p_lines$y4 <- p_lines$x+log10(2)

  main_fig<-
    ggplot(data = data,
         aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, group = runID))+
    theme_classic()+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 2, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 2, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 1, linewidth = 0.4, color = "grey30")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4, size = NULL, color = NULL, group = NULL, fill = NULL),
              linetype = 1, linewidth = 0.4, color = "grey30")+
    ggalt::geom_encircle(mapping = aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, fill = PFAA, group = PFAA),
                         alpha=0.3, s_shape=1, expand=0, size=1)+
    # stat_ellipse(mapping = aes(y = log_ngkg_re, x = Obs_logngkg, color = PFAA, fill = PFAA, group = PFAA))+
    # geom_smooth(method = "lm", se = FALSE, aes(color = NULL), color = "black", linewidth = 0.1)+
    # geom_abline(intercept = sum_est$mean_val[1], slope = sum_est$mean_val[2], size = 1, color = "black")+
    # geom_abline(intercept = sum_est$mean_val[1]+ sum_est$sd_val[1],
    #             slope = sum_est$mean_val[2]+ sum_est$sd_val[2],
    #             size = 0.2, linetype = 1, color = "red")+
    # geom_abline(intercept = sum_est$mean_val[1]- sum_est$sd_val[1],
    #             slope = sum_est$mean_val[2]- sum_est$sd_val[2],
    #             size = 0.2, linetype = 1, color = "red")+
    scale_color_manual(values = PFAA_cols)+
    scale_fill_manual(values = PFAA_cols)+
    xlab('Oberved log(ng/kg)')+
    ylab('Modeled log(ng/kg)')+
    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black", size = 1.5)+
    # theme(legend.position = "none")+
    geom_point(pch=19, size=1, alpha = 0.7)+
    ggtitle(plot_title)+
    ylim(0, 7)+
    xlim(0, 7)
  
  top_plots<-cowplot::plot_grid(p_cv, p_mean, pred_int, pred_b, nrow = 2)
  resid_plots<-cowplot::plot_grid( p_resid, p_resid1, ncol=1)
  full_plot<-cowplot::plot_grid(top_plots,
                                resid_plots, 
                                main_fig, nrow = 1)
  
  return(full_plot)
}

WG_observed_analysis<-rand_ecosystem_runs_analysis(data = as.data.frame(WG_observed[[1]]),
                             obs_metadata = metadata,
                             Site = "WG",
                             Season = "Fall",
                             plot_title = "WG - observed diet PFAA - 50 iter")
ggsave("./Figures/WG_observed_analysis.png", plot = WG_observed_analysis, width = 20, height = 6, units = "in")

WG_modeled_analysis<-rand_ecosystem_runs_analysis(data = as.data.frame(WG_modeled[[1]]),
                             obs_metadata = metadata,
                             Site = "WG",
                             Season = "Fall",
                             plot_title = "WG - modeled diet PFAA - 50 iter")
ggsave("./Figures/WG_modeled_analysis.png", plot = WG_modeled_analysis, width = 20, height = 6, units = "in")

# JBA 
JBAsp_observed_analysis<-rand_ecosystem_runs_analysis(data = as.data.frame(JBAsp_observed[[1]]),
                             obs_metadata = metadata,
                             Site = "JBA",
                             Season = "Spring",
                             plot_title = "JBA spring- observed diet PFAA - 50 iter")
ggsave("./Figures/JBAsp_observed_analysis.png", plot = JBAsp_observed_analysis, width = 20, height = 6, units = "in")


JBAsp_modeled_analysis<-rand_ecosystem_runs_analysis(data = as.data.frame(JBAsp_modeled[[1]]),
                             obs_metadata = metadata,
                             Site = "JBA",
                             Season = "Spring",
                             plot_title = "JBA spring - modeled diet PFAA - 50 iter")
ggsave("./Figures/JBAsp_modeled_analysis.png", plot = JBAsp_modeled_analysis, width = 20, height = 6, units = "in")


JBAsum_observed_analysis<-rand_ecosystem_runs_analysis(data = as.data.frame(JBAsum_observed[[1]]),
                             obs_metadata = metadata,
                             Site = "JBA",
                             Season = "Summer",
                             plot_title = "JBA summer- observed diet PFAA - 50 iter")
ggsave("./Figures/JBAsum_observed_analysis.png", plot = JBAsum_observed_analysis, width = 20, height = 6, units = "in")


JBAsum_modeled_analysis<-rand_ecosystem_runs_analysis(data = as.data.frame(JBAsum_modeled[[1]]),
                             obs_metadata = metadata,
                             Site = "JBA",
                             Season = "Summer",
                             plot_title = "JBA summer- modeled diet PFAA - 50 iter")
ggsave("./Figures/JBAsum_modeled_analysis.png", plot = JBAsum_modeled_analysis, width = 20, height = 6, units = "in")


# the line of identity is equivalent to:
# slope = 1 
# intercept = 0 
# Concordance Correlation Coefficient
# CCC(x, y, ci = "z-transform", conf.level = 0.95, na.rm = FALSE)
# intraclass correlation coefficient (ICC)

# Compare: 
# Same species different season



# sanity checks 
# ((AllData$Gill_uptake + AllData$Dietary_uptake) / AllData$`Total Elimination`) - AllData$C_B_re
# ((AllData$Gill_uptake + AllData$Dietary_uptake) / AllData$`Total Elimination`) - AllData$C_B_re

# ggplot(data = AllData,
#        aes(y = diff_Obs_mod_logngkg_re, x = Dietary_uptake,
#            color = PFAA, alpha = temperature, size = WB))+
#   geom_hline(yintercept = 0, linetype = "dashed")+
#   geom_point(pch=19)+
#   facet_wrap(PFAA~SppAlias, scale = "free")+
#   scale_color_manual(values = PFAA_cols)+
#   theme_classic()

# ggplot(data = AllData,
#        aes(y = diff_Obs_mod_logngkg_re, x = Gill_uptake,
#            color = PFAA, alpha = temperature, size = WB))+
#   geom_hline(yintercept = 0, linetype = "dashed")+
#   geom_point(pch=19)+
#   facet_wrap(PFAA~SppAlias, scale = "free")+
#   scale_color_manual(values = PFAA_cols)+
#   theme_classic()+
#   ggtitle("Bluegill")

# p3<-ggplot(data = AllData,
#        aes(x = Dietary_uptake, y = Diet, color = PFAA, size = WB,
#            shape = SppAlias, alpha = temperature, group = interaction(SppAlias, temperature)))+
#   geom_smooth(se=FALSE, linewidth = 0.3)+
#   geom_point()+
#   # geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
#   facet_wrap(PFAA~SppAlias, scales = "free")+
#   ylab(expression(PFAS~concentration~diet~(ng/g)))+
#   xlab(expression(Dietary~uptake~rate~(g[PFAS]/kg[fish]/d)))+
#   scale_color_manual(values = PFAA_cols)+
#   theme_bw()

# 
# plot3d(
#   z = AllData_blg$`temperature`, x = AllData_blg$`WB`, y = AllData_blg$`Dietary_uptake`,
#   type = 's',
#   col = AllData$`color`, 
#   zlab = "Temp", xlab="Mass", ylab="Dietary uptake")
# 
# plot3d(
#   z = AllData_blg$`temperature`, x = AllData_blg$`WB`, y = AllData_blg$`Gill_uptake`,
#   type = 's', size = 0.9,
#   col = AllData$`color`, 
#   zlab = "Temp", xlab="Mass", ylab="Gill uptake")
