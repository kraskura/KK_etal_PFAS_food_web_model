## run function:

# - spring vs summer 

library(dplyr)
library(tidyr)
library(ggplot2)
 
# source all functions to run all model runs
source("./Code/Sun et al reproduce/PFAS_model_Sun_etal_classes_KK.R") # Sets up classes
source("./Code/Sun et al reproduce/PFAS_model_Sun_etal_steadyStateMod_KK.R") # steady state model, needed for 'runEcosystemModel' it is called in Bioaccumulation model 
source("./Code/Sun et al reproduce/PFAS_model_Sun_etal_BioaccumulationMod.R") # steady state model, needed for 'runEcosystemModel'
source("./Code/runEcosystemModels.R") 
options(error = traceback)

# ***********************************************
# SECTION 1 all species, one site, one season : SUMMER --------------------------- 
#   site: JBA
#   species: loop through all species 
#   season: summer 
# ***********************************************

# 1. constants or original datasets; aply to all species loop-runs 
obs_BCF_JBA_summer_Allfish <- read.csv('./Data/JBA/BCF_BMF/all_species_BCF_JBA_summer.csv')
obs_BMF_JBA_summer_Allfish <- read.csv('./Data/JBA/BCF_BMF/all_species_BAF_JBA_summer.csv')
obs_size_JBA_summer_Allfish <- read.csv('./Data/JBA/BCF_BMF/all_species_massKG_JBA_summer.csv')
# guess_C_D_JBA_summer_Allfish <- read.csv('./Data/JBA_WG/all_species_C_D_JBA_summer_sedOnly.csv')
guess_C_D_JBA_summer_Allfish <- read.csv('./Data/JBA/BCF_BMF/all_species_C_D_JBA_summer.csv')

NgData = data.frame(kElimRatios = c(0.17, 0.74, 0.48, 0.74, 0.857, 0.95),
                    row.names = c('PFHxS', 'PFOS', 'PFOA', 'PFNA', 'PFDA', 'PFUA'))
Consoer_PFOSratio = 2.01/2.54 # clearance / total elimination ratio
Consoer_PFOAratio = 0.12/1.47 # clearance / total elimination ratio

settings_BCF = new_Settings(
    chooseStudyType = 'Martin BCF',
    chooseDiet= 'zero',
    chooseEd= 'Martin',
    chooseRenal= 'off'
)
settings_BMF = new_Settings(
    chooseStudyType= 'Martin BMF',
    chooseDiet= 'Martin BMF',
    chooseEd= 'Martin',
    chooseRenal= 'off'
)
settings_BMFd = new_Settings(
    chooseStudyType= 'Martin BMF',
    chooseDiet= 'default',
    chooseEd= 'Martin',
    chooseRenal= 'off'
)

# !! change size in organismal data for every loop-run  
BCF_inputFiles_JBA_summer = list(
    'numSpecies' = 2,
    'oceanData' = read.csv('./Data/JBA/BCF_BMF/oceanData_BCFBMF_JBA_summer.csv', row.names = "oceanParameter"),
    'chemicalParams'= read.csv('./Data/JBA/BCF_BMF/chemicalParameters_BCFBMF_JBA.csv', row.names = "chemicalParameter"),
    # 'chemicalData'= read.csv('./Data/JBA/BCF_BMF/chemicalData_BCF_JBA_summer_max.csv', row.names = "chemicalParameter"),
    'chemicalData'= read.csv('./Data/JBA/BCF_BMF/chemicalData_BCF_JBA_summer_med.csv', row.names = "chemicalParameter"),
    'organismData'= read.csv('./Data/JBA/BCF_BMF/organismData_BMFBCF_JBA.csv', row.names = "X"),
    'foodWebData'= read.csv('./Data/JBA/BCF_BMF/foodWebTable_BCFBMF_JBA_fish.csv', row.names = "X"))
BMF_inputFiles_JBA_summer = list(
    'numSpecies'= 2,
    'oceanData'= read.csv('./Data/JBA/BCF_BMF/oceanData_BCFBMF_JBA_summer.csv', row.names = "oceanParameter"),
    'chemicalParams'= read.csv('./Data/JBA/BCF_BMF/chemicalParameters_BCFBMF_JBA.csv', row.names = "chemicalParameter"),
    # 'chemicalData'= read.csv('./Data/JBA/BCF_BMF/chemicalData_BMF_JBA_summer_max.csv', row.names = "chemicalParameter"),
    'chemicalData'= read.csv('./Data/JBA/BCF_BMF/chemicalData_BMF_JBA_summer_med.csv', row.names = "chemicalParameter"),
    'organismData'= read.csv('./Data/JBA/BCF_BMF/organismData_BMFBCF_JBA.csv', row.names = "X"),
    'foodWebData'= read.csv('./Data/JBA/BCF_BMF/foodWebTable_BCFBMF_JBA_fish.csv', row.names = "X"))

PFAA_BCFs_ALLspecies<-NULL
PFAA_BMFs_ALLspecies<-NULL

for (i in 1:nrow(obs_BCF_JBA_summer_Allfish)){
    
    species <- obs_BCF_JBA_summer_Allfish$Species[i]
    print(species)
    BCF_inputFiles_JBA_summer$organismData[1,"Fish"] <- obs_size_JBA_summer_Allfish$mass_kg[i]
    BMF_inputFiles_JBA_summer$chemicalData["C_D",] <- guess_C_D_JBA_summer_Allfish[i,c(2:8)]

    observedData_BCF <- data.frame(logBCF_obs_WBcalc = c(log10(obs_BCF_JBA_summer_Allfish$PFHxS[i]),
                                                         log10(obs_BCF_JBA_summer_Allfish$PFOS[i]),
                                                         log10(obs_BCF_JBA_summer_Allfish$PFOA[i]),
                                                         log10(obs_BCF_JBA_summer_Allfish$PFNA[i]),
                                                         log10(obs_BCF_JBA_summer_Allfish$PFDA[i]),
                                                         log10(obs_BCF_JBA_summer_Allfish$PFUA[i])),
                                   BCF_obs_carcass = c(obs_BCF_JBA_summer_Allfish$PFHxS[i],
                                                         obs_BCF_JBA_summer_Allfish$PFOS[i],
                                                         obs_BCF_JBA_summer_Allfish$PFOA[i],
                                                         obs_BCF_JBA_summer_Allfish$PFNA[i],
                                                         obs_BCF_JBA_summer_Allfish$PFDA[i],
                                                         obs_BCF_JBA_summer_Allfish$PFUA[i]), # Large mouth bass, whole body here. 
                                   row.names = c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA"))
     print(observedData_BCF$BCF_obs_carcass)
    if(all(is.nan(unlist(observedData_BCF[1,])))){
        next
    }
    
    observedData_BMF <- data.frame(BAF_obs = c(obs_BMF_JBA_summer_Allfish$PFHxS[i],
                                                 obs_BMF_JBA_summer_Allfish$PFOS[i],
                                                 obs_BMF_JBA_summer_Allfish$PFOA[i],
                                                 obs_BMF_JBA_summer_Allfish$PFNA[i],
                                                 obs_BMF_JBA_summer_Allfish$PFDA[i],
                                                 obs_BMF_JBA_summer_Allfish$PFUA[i]),	# Large mouth bass				, 
                                   row.names = c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA") ) # reported as BAF, not BMF
    # Get BCF and BMF data 
    # PFAA_BMFs_out <- runEcosystemModel(settings = settings_BMF,
    #                   inputFiles_list = BMF_inputFiles_JBA_summer,
    #                   PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
    #                   BMF_calc = TRUE,
    #                   NgData = NgData,
    #                   Consoer_PFOSratio = Consoer_PFOSratio, 
    #                   Consoer_PFOAratio = Consoer_PFOAratio,
    #                   observedData_BCF_BMF = observedData_BMF)
    
    PFAA_BCFs_out <- runEcosystemModel(settings = settings_BCF,
                      inputFiles_list = BCF_inputFiles_JBA_summer,
                      PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                      BCF_calc = TRUE,
                      NgData = NgData,
                      Consoer_PFOSratio = Consoer_PFOSratio, 
                      Consoer_PFOAratio = Consoer_PFOAratio,
                      observedData_BCF_BMF = observedData_BCF)

    PFAA_BCFs<-PFAA_BCFs_out[[1]]
    # model_BCF<-PFAA_BCFs_out[[2]]
    # df_BCFreg<-PFAA_BCFs_out[[3]]
    # pred_BCF<-PFAA_BCFs_out[[4]]
    PFAA_BCFs$BCF_BMF<-"BCF"
    PFAA_BCFs$Species<-species
    PFAA_BCFs$season<-"summer"
    PFAA_BCFs$W_B<-obs_size_JBA_summer_Allfish$mass_kg[i]
    
    # PFAA_BMFs<-PFAA_BMFs_out[[1]]
    # # model_BMF<-PFAA_BMFs_out[[2]]
    # # df_BMFreg<-PFAA_BMFs_out[[3]]
    # # pred_BMF<-PFAA_BMFs_out[[4]]
    # PFAA_BMFs$BCF_BMF<-"BMF"
    # PFAA_BMFs$Species<-species
    # PFAA_BMFs$W_B<-obs_size_JBA_summer_Allfish$mass_kg[i]
    # PFAA_BMFs$season<-"summer"

    if(i == 1){
        PFAA_BCFs_ALLspecies<-PFAA_BCFs
        # PFAA_BMFs_ALLspecies<-PFAA_BMFs
    }else{
        PFAA_BCFs_ALLspecies<-rbind(PFAA_BCFs_ALLspecies, PFAA_BCFs)
        # PFAA_BMFs_ALLspecies<-rbind(PFAA_BMFs_ALLspecies, PFAA_BMFs)
    }
}

PFAA_BCFs_ALLspecies_summer <- PFAA_BCFs_ALLspecies
# PFAA_BMFs_ALLspecies_summer <- PFAA_BMFs_ALLspecies
 

# ***********************************************
# SECTION 2 all species, one site, one season : SPRING --------------------------- 
#   site: JBA
#   species: loop through all species 
#   season: spring
# ***********************************************
 
# 1. constants or original data sets; apply to all species loop-runs 
obs_BCF_JBA_spring_Allfish <- read.csv('./Data/JBA/BCF_BMF/all_species_BCF_JBA_spring.csv')
obs_BMF_JBA_spring_Allfish <- read.csv('./Data/JBA/BCF_BMF/all_species_BAF_JBA_spring.csv')
obs_size_JBA_spring_Allfish <- read.csv('./Data/JBA/BCF_BMF/all_species_massKG_JBA_spring.csv')
# guess_C_D_JBA_spring_Allfish <- read.csv('./Data/JBA_WG/all_species_C_D_JBA_spring_sedOnly.csv')
guess_C_D_JBA_spring_Allfish <- read.csv('./Data/JBA/BCF_BMF/all_species_C_D_JBA_spring.csv')


NgData = data.frame(kElimRatios = c(0.17, 0.74, 0.48, 0.74, 0.857, 0.95),
                    row.names = c('PFHxS', 'PFOS', 'PFOA', 'PFNA', 'PFDA', 'PFUA'))
Consoer_PFOSratio = 2.01/2.54 # clearance / total elimination ratio
Consoer_PFOAratio = 0.12/1.47 # clearance / total elimination ratio

settings_BCF = new_Settings(
    chooseStudyType = 'Martin BCF',
    chooseDiet= 'zero',
    chooseEd= 'Martin',
    chooseRenal= 'off'
)
settings_BMF = new_Settings(
    chooseStudyType= 'Martin BMF',
    chooseDiet= 'Martin BMF',
    chooseEd= 'Martin',
    chooseRenal= 'off'
)

# !! change size in organismal data for every loop-run  
BCF_inputFiles_JBA_spring = list(
    'numSpecies' = 2,
    'oceanData' = read.csv('./Data/JBA/BCF_BMF/oceanData_BCFBMF_JBA_spring.csv', row.names = "oceanParameter"),
    'chemicalParams'= read.csv('./Data/JBA/BCF_BMF/chemicalParameters_BCFBMF_JBA.csv', row.names = "chemicalParameter"),
    # 'chemicalData'= read.csv('./Data/JBA/BCF_BMF/chemicalData_BCF_JBA_spring_max.csv', row.names = "chemicalParameter"),
    'chemicalData'= read.csv('./Data/JBA/BCF_BMF/chemicalData_BCF_JBA_spring_med.csv', row.names = "chemicalParameter"),
    'organismData'= read.csv('./Data/JBA/BCF_BMF/organismData_BMFBCF_JBA.csv', row.names = "X"),
    'foodWebData'= read.csv('./Data/JBA/BCF_BMF/foodWebTable_BCFBMF_JBA_fish.csv', row.names = "X"))
BMF_inputFiles_JBA_spring = list(
    'numSpecies'= 2,
    'oceanData'= read.csv('./Data/JBA/BCF_BMF/oceanData_BCFBMF_JBA_spring.csv', row.names = "oceanParameter"),
    'chemicalParams'= read.csv('./Data/JBA/BCF_BMF/chemicalParameters_BCFBMF_JBA.csv', row.names = "chemicalParameter"),
    'chemicalData'= read.csv('./Data/JBA/BCF_BMF/chemicalData_BMF_JBA_spring_med.csv', row.names = "chemicalParameter"),
    # 'chemicalData'= read.csv('./Data/JBA/BCF_BMF/chemicalData_BMF_JBA_spring_max.csv', row.names = "chemicalParameter"),
    'organismData'= read.csv('./Data/JBA/BCF_BMF/organismData_BMFBCF_JBA.csv', row.names = "X"),
    'foodWebData'= read.csv('./Data/JBA/BCF_BMF/foodWebTable_BCFBMF_JBA_fish.csv', row.names = "X"))

PFAA_BCFs_ALLspecies<-NULL
PFAA_BMFs_ALLspecies<-NULL

for (i in 1:nrow(obs_BCF_JBA_spring_Allfish)){
    
    species <- obs_BCF_JBA_spring_Allfish$Species[i]
    print(species)
    BCF_inputFiles_JBA_spring$organismData[1,"Fish"] <- obs_size_JBA_spring_Allfish$mass_kg[i]
    BMF_inputFiles_JBA_spring$chemicalData["C_D",] <- guess_C_D_JBA_spring_Allfish[i,c(2:8)]

    observedData_BCF <- data.frame(logBCF_obs_WBcalc = c(log10(obs_BCF_JBA_spring_Allfish$PFHxS[i]),
                                                         log10(obs_BCF_JBA_spring_Allfish$PFOS[i]),
                                                         log10(obs_BCF_JBA_spring_Allfish$PFOA[i]),
                                                         log10(obs_BCF_JBA_spring_Allfish$PFNA[i]),
                                                         log10(obs_BCF_JBA_spring_Allfish$PFDA[i]),
                                                         log10(obs_BCF_JBA_spring_Allfish$PFUA[i])),
                                   BCF_obs_carcass = c(obs_BCF_JBA_spring_Allfish$PFHxS[i],
                                                         obs_BCF_JBA_spring_Allfish$PFOS[i],
                                                         obs_BCF_JBA_spring_Allfish$PFOA[i],
                                                         obs_BCF_JBA_spring_Allfish$PFNA[i],
                                                         obs_BCF_JBA_spring_Allfish$PFDA[i],
                                                         obs_BCF_JBA_spring_Allfish$PFUA[i]), # Large mouth bass, whole body here. 
                                   row.names = c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA"))
     print(observedData_BCF$BCF_obs_carcass)
    if(all(is.nan(unlist(observedData_BCF[1,])))){
        next
    }
    
    observedData_BMF <- data.frame(BAF_obs = c(obs_BMF_JBA_spring_Allfish$PFHxS[i],
                                                 obs_BMF_JBA_spring_Allfish$PFOS[i],
                                                 obs_BMF_JBA_spring_Allfish$PFOA[i],
                                                 obs_BMF_JBA_spring_Allfish$PFNA[i],
                                                 obs_BMF_JBA_spring_Allfish$PFDA[i],
                                                 obs_BMF_JBA_spring_Allfish$PFUA[i]),	# Large mouth bass				, 
                                   row.names = c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA") ) # reported as BAF, not BMF
    # Get BCF and BMF data 
    PFAA_BMFs_out <- runEcosystemModel(settings = settings_BMF,
                      inputFiles_list = BMF_inputFiles_JBA_spring,
                      PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                      BMF_calc = TRUE,
                      NgData = NgData,
                      Consoer_PFOSratio = Consoer_PFOSratio, 
                      Consoer_PFOAratio = Consoer_PFOAratio,
                      observedData_BCF_BMF = observedData_BMF)
    
    PFAA_BCFs_out <- runEcosystemModel(settings = settings_BCF,
                      inputFiles_list = BCF_inputFiles_JBA_spring,
                      PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                      BCF_calc = TRUE,
                      NgData = NgData,
                      Consoer_PFOSratio = Consoer_PFOSratio, 
                      Consoer_PFOAratio = Consoer_PFOAratio,
                      observedData_BCF_BMF = observedData_BCF)

    PFAA_BCFs<-PFAA_BCFs_out[[1]]
    # model_BCF<-PFAA_BCFs_out[[2]]
    # df_BCFreg<-PFAA_BCFs_out[[3]]
    # pred_BCF<-PFAA_BCFs_out[[4]]
    PFAA_BCFs$BCF_BMF<-"BCF"
    PFAA_BCFs$Species<-species
    PFAA_BCFs$season<-"spring"
    PFAA_BCFs$W_B<-obs_size_JBA_summer_Allfish$mass_kg[i]
    
    PFAA_BMFs<-PFAA_BMFs_out[[1]]
    # model_BMF<-PFAA_BMFs_out[[2]]
    # df_BMFreg<-PFAA_BMFs_out[[3]]
    # pred_BMF<-PFAA_BMFs_out[[4]]
    PFAA_BMFs$BCF_BMF<-"BMF"
    PFAA_BMFs$Species<-species
    PFAA_BMFs$W_B<-obs_size_JBA_summer_Allfish$mass_kg[i]
    PFAA_BMFs$season<-"spring"

    if(i == 1){
        PFAA_BCFs_ALLspecies<-PFAA_BCFs
        PFAA_BMFs_ALLspecies<-PFAA_BMFs
    }else{
        PFAA_BCFs_ALLspecies<-rbind(PFAA_BCFs_ALLspecies, PFAA_BCFs)
        PFAA_BMFs_ALLspecies<-rbind(PFAA_BMFs_ALLspecies, PFAA_BMFs)
    }
}

PFAA_BCFs_ALLspecies<-rbind(PFAA_BCFs_ALLspecies, PFAA_BCFs_ALLspecies_summer)
PFAA_BMFs_ALLspecies<-rbind(PFAA_BMFs_ALLspecies, PFAA_BMFs_ALLspecies_summer)
colnames(PFAA_BMFs_ALLspecies)
colnames(PFAA_BMFs_ALLspecies_summer)

# Plotting -------
 
# plot_REratios <- PFAA_BCFs[,'kr/totelim', drop =FALSE] # rounded to 3 decimal places BCF and BMF are the same
# plot_REratios$PFCA_CNo_est <- c(NaN, NaN, NaN, 8, 10, 11)
# plot_REratios$PFSA_CNo_est <- c(6, NaN, NaN, NaN, NaN, NaN)
# plot_REratios$PFCA_CNo <- c(NaN, NaN, 7, NaN, NaN, NaN)
# plot_REratios$PFSA_CNo <- c(NaN, 8, NaN, NaN, NaN, NaN)

rbind()
PFAA_cols<-c('PFHxS' = 'blueviolet',
             'PFOS'  = 'cornflowerblue',
             'PFOA'  = 'orange',
             'PFDA'  = 'deeppink',
             'PFUA'  = 'limegreen', 
             'PFDoDAc' = "grey40")

p_lines <- data.frame('x' = c(-3, 7),
                      'y' = c(-3, 7))
p_lines$y1 <- p_lines$x-1
p_lines$y2 <- p_lines$x+1
p_lines$y3 <- p_lines$x-log10(2)
p_lines$y4 <- p_lines$x+log10(2)

PFAA_BCFs$logBCF_obs_WBcalc-PFAA_BCFs$logBCF_obs_carcass

p1_noRen_Elim_BCF<-
  ggplot(data = PFAA_BCFs_ALLspecies, aes(x = logBCF, y = logBCF_obs_carcass, fill = PFAA))+
  geom_line(data = p_lines, mapping = aes(x = x, y = y),
            linetype = 1, linewidth = 1, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y1),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y2),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y3),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y4),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_point(pch=21, size = 3)+
  # geom_point(data=PFAA_BMFs_ALLspecies, aes(x = logBMF, y = logBAF_obs,
  #                                           fill = PFAA, size = W_B), pch = 23)+
  facet_grid(.~season)+
  theme_classic()+
  scale_fill_manual(values = PFAA_cols)+
  scale_color_manual(values = PFAA_cols)+
  xlab('Modeled log BCF (L/kg)')+
  ylab('Observed log BCF (L/kg)')+
  coord_cartesian(ylim=c(-2, 4), xlim = c(-2,4))+
  theme(legend.position = "none")+
  ggtitle("Without renal elimination")

p1_noRen_Elim_BMF<-
  ggplot(data=PFAA_BMFs_ALLspecies, aes(x = logBMF, y = logBAF_obs, fill = PFAA))+
  geom_line(data = p_lines, mapping = aes(x = x, y = y),
            linetype = 1, linewidth = 1, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y1),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y2),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y3),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y4),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_point(pch=21, size = 3)+
  # geom_point(, pch = 23)+
  facet_grid(.~season)+
  theme_classic()+
  scale_fill_manual(values = PFAA_cols)+
  scale_color_manual(values = PFAA_cols)+
  xlab('Modeled log BMF (kg/kg)')+
  ylab('Observed log BMF (kg/kg)')+
  coord_cartesian(ylim=c(-2, 4), xlim = c(-2,4))+
  theme(legend.position = "none")+
  ggtitle("Without renal elimination")

p2_wRen_Elim_BCF<-ggplot(data=PFAA_BCFs_ALLspecies, aes(x = logBCF_recalc, y = logBCF_obs_carcass, fill = PFAA))+
  geom_line(data = p_lines, mapping = aes(x = x, y = y),
            linetype = 1, linewidth = 1, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y1),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y2),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y3),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y4),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_point(pch=21, size = 3)+
  facet_grid(.~season)+
  # geom_point(data=PFAA_BMFs, aes(x = logBMF_recalc, y = logBAF_obs, fill = PFAA, size = W_B), pch=23)+
  scale_fill_manual(values = PFAA_cols)+
  scale_color_manual(values = PFAA_cols)+
  xlab('Modeled log BCF (L/kg)')+
  ylab('Observed log BCF (L/kg)')+
  coord_cartesian(ylim=c(-2, 4), xlim = c(-2,4))+
  theme_classic()+
  ggtitle("With renal elimination")+
  theme(legend.position = "none")


p2_wRen_Elim_BMF<-ggplot(data=PFAA_BMFs_ALLspecies, aes(x = logBMF_recalc, y = logBAF_obs, fill = PFAA))+
  geom_line(data = p_lines, mapping = aes(x = x, y = y),
            linetype = 1, linewidth = 1, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y1),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y2),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y3),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_line(data = p_lines, mapping = aes(x = x, y = y4),
            linetype = 3, linewidth = 0.6, inherit.aes = FALSE)+
  geom_point(pch=21, size = 3)+
  facet_grid(.~season)+
  scale_fill_manual(values = PFAA_cols)+
  scale_color_manual(values = PFAA_cols)+
  xlab('Modeled log BMF (kg/kg)')+
  ylab('Observed log BMF (kg/kg)')+
  coord_cartesian(ylim=c(-2, 4), xlim = c(-2,4))+
  theme_classic()+
  ggtitle("With renal elimination")+
  theme(legend.position = "none")



cowplot::plot_grid(p1_noRen_Elim_BCF,p2_wRen_Elim_BCF,
                   p1_noRen_Elim_BMF, p2_wRen_Elim_BMF,
                   nrow = 2,
                   ncol =2,
                   labels = "AUTO",
                   rel_widths = c(1, 1))


p_bar_BCF<-  ggplot(data = PFAA_BCFs_ALLspecies, aes(x = PFAA, y = logBCF_obs_carcass,
                                         fill = PFAA, color = PFAA))+ # observed 
  geom_bar(stat = "identity", fill = "white", width = 0.1, position = position_nudge(x=0.2))+
  geom_bar(aes(x = PFAA, y = logBCF_recalc, fill = PFAA), stat = "identity", alpha = 0.5, position = position_nudge(x=0.1), width = 0.1)+ # modeled w/ renal
  geom_bar(aes(x = PFAA, y = logBCF, fill = PFAA), stat = "identity", alpha = 0.9, position = position_nudge(x=0), width = 0.1)+ # modeled wo/ Renal 
  theme_classic()+
  coord_flip()+
  facet_grid(season~Species)+
  scale_shape_manual(values = c("BCF" = 21, "BMF" = 23))+
  scale_fill_manual(values = PFAA_cols)+
  scale_color_manual(values = PFAA_cols)+
  ylab('log BCF (L/kg)')+
  xlab('')+
  theme(legend.position = "none")
 
p_bar_BMF<- ggplot(data = PFAA_BMFs_ALLspecies, aes(x = PFAA, y = logBAF_obs,
                                         fill = PFAA, color = PFAA))+ # observed 
  geom_bar(stat = "identity", fill = "white", width = 0.1, position = position_nudge(x=0.2))+
  geom_bar(aes(x = PFAA, y = logBMF_recalc, fill = PFAA), stat = "identity", alpha = 0.5, position = position_nudge(x=0.1), width = 0.1)+ # modeled w/ renal
  geom_bar(aes(x = PFAA, y = logBMF, fill = PFAA), stat = "identity", alpha = 0.9, position = position_nudge(x=0), width = 0.1)+ # modeled wo/ Renal 
  theme_classic()+
  coord_flip()+
  facet_grid(season~Species)+
  scale_shape_manual(values = c("BCF" = 21, "BMF" = 23))+
  scale_fill_manual(values = PFAA_cols)+
  scale_color_manual(values = PFAA_cols)+
  ylab('log BMF (kg/kg)')+
  xlab('')+
  theme(legend.position = "none")

colnames(PFAA_BCFs_ALLspecies)<-paste(colnames(PFAA_BCFs_ALLspecies),"BCF",sep="_")
colnames(PFAA_BMFs_ALLspecies)<-paste(colnames(PFAA_BMFs_ALLspecies),"BMF",sep="_")

data<-cbind(PFAA_BCFs_ALLspecies, PFAA_BMFs_ALLspecies)
# 
# ggplot(data, aes(MB_BCF, MB_BMF, fill = PFAA_BMF, shape = season_BF))+
#     geom_point(pch=21, size = 3)+
#     xlim(0, 5)+
#     ylim(0,5)
# 
# 
# 









