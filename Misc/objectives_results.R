## ----setup1, include=FALSE------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE,
                      message = FALSE,
                      warning = FALSE,
                      fig.width=8.5,
                      fig.height=4)


## ----source, warning=FALSE, message=FALSE, echo=FALSE---------------------------------------------------

library(kableExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)
library(tibble) # needed to get named list of dataframes (in data_tables.R)
library(readxl)
library(ggformat2)
library(here)
here::i_am(path = "./Code/R results/objectives_results.Rmd")
# here()

# functions to run food web model
source(here("Code", "R model", "PFAS_bioaccum_v2.R"))
source(here("Code", "R model", "PFAS_classes_v2.R"))
source(here("Code", "R model", "PFAS_steadyState_v2.R"))
source(here("Code", "R model", "runEcosystemModels.R")) # custom written
source(here("Code", "R model", "runErrorEstimates.R"))
source(here("Code", "R model", "data_tables.R"))


# source("../Code/PFAS_tissuePart_v2.R")

# source data tables for each objective
kRTable = read.csv(here("Data", 'export_krTable.csv'), row.names = "X", header = T)
colnames(kRTable)<-c("kr/kb", "PFAA", "chemID")

PFAS_List <- c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')
PFAS_labels <- c('C6 PFSA','C8 PFSA','C8 PFCA', 'C9 PFCA','C10 PFCA','C11 PFCA')
PFAA_cols <- c('PFHxS' = '#332288', # 
             'PFOS' = '#882255', # 
             'PFOA' = '#117733',
             'PFNA' = '#88CCEE',
             'PFDA' = 'orange',
             'PFUA' = '#AA4499')

p_lines <- data.frame('x' = c(0, 7),
                      'y' = c(0, 7))
p_lines$y1 <- p_lines$x-1
p_lines$y2 <- p_lines$x+1
p_lines$y3 <- p_lines$x-log10(2)
p_lines$y4 <- p_lines$x+log10(2)



## ----willow grove all data------------------------------------------------------------------------------
source(here("Code", "R data tables", "WG.R"))

# PFAS are known in each diet item
settings_obsDiet = new_Settings(
  chooseDiet = 'forced',
  chooseKoc = 'Koc',
  chooseDmw = "Droge",
  chooseEw = "empirical"
)

# used below in sed:water est part 
WG_inputFiles_list<-inputFiles_list
water_WG<-d.water %>%  # get the medians from the data
  select(PFHxS.m, PFOS.m, PFOA.m, PFNA.m) %>% 
  as.data.frame()
water_WG<-water_WG/1000
colnames(water_WG)<-c("PFHxS", "PFOS", "PFOA", "PFNA") 
sed_WG<-d.sed %>%  # get the medians from the data
  select(PFHxS.m, PFOS.m, PFOA.m, PFNA.m) %>% 
  as.data.frame()
colnames(sed_WG)<-c("PFHxS", "PFOS", "PFOA", "PFNA") 
obs_ratio_WG<-sed_WG/water_WG


# run the bioaccumulation model: 
AllDataWG<-runEcosystemModel(settings = settings_obsDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR',
                                    "Gill_uptake","Dietary_uptake", "Sediment_uptake",
                                    "Gill_up%", "Diet_up%", 'Sed_up%',
                                    "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = inputFiles_list$max_dietData)

# PFAS in each diet item are calculated from trophic transfer
settings_modelDiet = new_Settings(
  chooseDiet = 'default',
  chooseKoc = 'Koc',
  chooseDmw = "Droge", # default = "Droge", or use "calculated" 
  chooseEw = "empirical", # default = 'empirical', or use "beta#", where:
  # must be in format: beta#, #beta [# = numeric var, no spaces, no special char]"
)

AllDataWG_fw<-runEcosystemModel(settings = settings_modelDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR',
                                    "Gill_uptake","Dietary_uptake", "Sediment_uptake",
                                    "Gill_up%", "Diet_up%", 'Sed_up%',
                                    "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = inputFiles_list$max_dietData)

# PFAS in each diet item are calculated from trophic transfer
settings_modelDiet = new_Settings(
  chooseDiet = 'zero',
  chooseKoc = 'Koc',
  chooseDmw = "Droge", # default = "Droge", or use "calculated" 
  chooseEw = "empirical", # default = 'empirical', or use "beta#", where:
  # must be in format: beta#, #beta [# = numeric var, no spaces, no special char]"
)

AllDataWG_zero<-runEcosystemModel(settings = settings_modelDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR',
                                    "Gill_uptake","Dietary_uptake", "Sediment_uptake",
                                    "Gill_up%", "Diet_up%", 'Sed_up%',
                                    "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = inputFiles_list$max_dietData)



## ----willow grove same date samping---------------------------------------------------------------------
source(here("Code", "R data tables", "WG_sameDate.R"))

# PFAS are known in each diet item
settings_obsDiet = new_Settings(
  chooseDiet = 'forced',
  chooseKoc = 'Koc',
  chooseDmw = "Droge",
  chooseEw = "empirical"
)

# run the bioaccumulation model: 
AllDataWG2<-runEcosystemModel(settings = settings_obsDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR', "Gill_uptake",
                                    "Gill_up%", "Diet_up%",
                                    "Dietary_uptake", "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = inputFiles_list$max_dietData)

# PFAS in each diet item are calculated from trophic transfer
settings_modelDiet = new_Settings(
  chooseDiet = 'default',
  chooseKoc = 'Koc',
  chooseDmw = "Droge", # default = "Droge", or use "calculated" 
  chooseEw = "empirical", # default = 'empirical', or use "beta#", where:
  # must be in format: beta#, #beta [# = numeric var, no spaces, no special char]"
)

AllDataWG_fw2<-runEcosystemModel(settings = settings_modelDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR', "Gill_uptake",
                                    "Gill_up%", "Diet_up%",
                                    "Dietary_uptake", "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = inputFiles_list$max_dietData)



## ----wg model bias tables, echo = FALSE, message=FALSE, warning=FALSE-----------------------------------

wg1_error <- estimate_error(MetricData = AllDataWG, DataID = "wg1")
wg1_error_fw <- estimate_error(MetricData = AllDataWG_fw, DataID = "wg1")

wg1_error2 <- estimate_error(MetricData = AllDataWG2, DataID = "wg1")
wg1_error_fw2 <- estimate_error(MetricData = AllDataWG_fw2, DataID = "wg1")
# MBj, MBi, MB, rsq,  r2_pfaa, r2_species


## ----willow grove data output and plots, fig.cap= "Figure. 1"-------------------------------------------
AllDataWG <- AllDataWG[!c(AllDataWG$PFAA == "PFUA" | AllDataWG$PFAA == "PFDA"),]
AllDataWG_fw <- AllDataWG_fw[!c(AllDataWG_fw$PFAA == "PFUA" | AllDataWG_fw$PFAA == "PFDA"),]
AllDataWG_zero <- AllDataWG_zero[!c(AllDataWG_zero$PFAA == "PFUA" | AllDataWG_zero$PFAA == "PFDA"),]
AllDataWG2 <- AllDataWG[!c(AllDataWG$PFAA == "PFUA" | AllDataWG$PFAA == "PFDA"),]
AllDataWG_fw2 <- AllDataWG_fw[!c(AllDataWG_fw$PFAA == "PFUA" | AllDataWG_fw$PFAA == "PFDA"),]

p1<-ggplot(data = AllDataWG, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo, ymax = Obs_logngkg + Obs_logngkg_Up), linewidth = 0.5) +
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo, xmax = log_ngkg_re + logngkg_re_Up), linewidth = 0.5) +
    geom_point(pch = 19, size = 3)+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Known PFAS in diet')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Measured log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(2,6.5), xlim = c(2,6.5))+
    theme(legend.position = "none")
# p1_fish

p2<-ggplot(data = AllDataWG_fw, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo, ymax = Obs_logngkg + Obs_logngkg_Up), linewidth = 0.5) +
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo, xmax = log_ngkg_re + logngkg_re_Up), linewidth = 0.5) +
    geom_point(pch = 19, size = 3)+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Modeled trophic transfer of PFAA')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Measured log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(0,6.5), xlim = c(0,6.5))+
    # theme(legend.position = "right")
    theme(legend.position = "none")


p3<-ggplot(data = AllDataWG_zero, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo, ymax = Obs_logngkg + Obs_logngkg_Up), linewidth = 0.5) +
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo, xmax = log_ngkg_re + logngkg_re_Up), linewidth = 0.5) +
    geom_point(pch = 19, size = 3)+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'ZERO diet')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Measured log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(2,6.5), xlim = c(2,6.5))+
    theme(legend.position = "none")


p1.2<-ggplot(data = AllDataWG2, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo, ymax = Obs_logngkg + Obs_logngkg_Up), linewidth = 0.5) +
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo, xmax = log_ngkg_re + logngkg_re_Up), linewidth = 0.5) +
    geom_point(pch = 19, size = 3)+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Known PFAS in diet: SAME DAY samples')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Measured log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(2,6.5), xlim = c(2,6.5))+
    theme(legend.position = "none")
# p1_fish

p2.2<-ggplot(data = AllDataWG_fw2, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo, ymax = Obs_logngkg + Obs_logngkg_Up), linewidth = 0.5) +
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo, xmax = log_ngkg_re + logngkg_re_Up), linewidth = 0.5) +
    geom_point(pch = 19, size = 3)+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Modeled trophic transfer of PFAA: SAME DAY samples')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Measured log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(0,6.5), xlim = c(0,6.5))+
    # theme(legend.position = "right")
    theme(legend.position = "none")

# p1_fish_re
legend_plot <- ggplot(data = AllDataWG_fw,
         aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
  geom_point()+
  scale_color_manual(values = PFAA_cols)+
  theme(legend.title = element_blank())+
  theme_bw()

legend<-cowplot::get_legend(legend_plot)
fig1<-cowplot::plot_grid(p1, p2, p3,
                         p1.2, p2.2, legend,
                   ncol = 3,
                   nrow = 2,
                   rel_widths = c(1,1, 1),
        labels = c("A", "B", "C", "D", "E", ""))
fig1



## ----jba spring all data--------------------------------------------------------------------------------
source(here("Code", "R data tables", "JBA_spring.R"))

settings_obsDiet = new_Settings(
  chooseDiet = 'forced',
)

# used below for sed:weater estimates
JBAsp_inputFiles_list<-inputFiles_list
water_JBAsp<-d.water %>%  # get the medians from the data
  select(PFHxS.m, PFOS.m, PFOA.m, PFNA.m) %>% 
  as.data.frame()
water_JBAsp<-water_JBAsp/1000
colnames(water_JBAsp)<-c("PFHxS", "PFOS", "PFOA", "PFNA") 
sed_JBAsp<-d.sed %>%  # get the medians from the data
  select(PFHxS.m, PFOS.m, PFOA.m, PFNA.m) %>% 
  as.data.frame()
colnames(sed_JBAsp)<-c("PFHxS", "PFOS", "PFOA", "PFNA") 
obs_ratio_JBAsp<-sed_JBAsp/water_JBAsp


AllDataJBAsp<-runEcosystemModel(settings = settings_obsDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR',
                                    "Gill_uptake","Dietary_uptake", "Sediment_uptake",
                                    "Gill_up%", "Diet_up%", 'Sed_up%',
                                    "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = inputFiles_list$max_dietData)


settings_modelDiet = new_Settings(
  chooseDiet = 'default'
)

AllDataJBAsp_fw<-runEcosystemModel(settings = settings_modelDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR',
                                    "Gill_uptake","Dietary_uptake", "Sediment_uptake",
                                    "Gill_up%", "Diet_up%", 'Sed_up%',
                                    "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = inputFiles_list$max_dietData)

AllDataJBAsp <- AllDataJBAsp[!c(AllDataJBAsp$PFAA == "PFUA" ),]
AllDataJBAsp_fw <- AllDataJBAsp_fw[!c(AllDataJBAsp_fw$PFAA == "PFUA" ),]




## ----jba spring same day data---------------------------------------------------------------------------
source(here("Code", "R data tables", "JBA_spring_sameDate.R"))

settings_obsDiet = new_Settings(
  chooseDiet = 'forced',
)

AllDataJBAsp2<-runEcosystemModel(settings = settings_obsDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR', "Gill_uptake",
                                    "Gill_up%", "Diet_up%",
                                    "Dietary_uptake", "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = inputFiles_list$max_dietData)


settings_modelDiet = new_Settings(
  chooseDiet = 'default'
)

AllDataJBAsp_fw2<-runEcosystemModel(settings = settings_modelDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR', "Gill_uptake",
                                    "Gill_up%", "Diet_up%", "Dietary_uptake",
                                    "Total Elimination", "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = inputFiles_list$max_dietData)

AllDataJBAsp2 <- AllDataJBAsp2[!c(AllDataJBAsp2$PFAA == "PFUA" ),]
AllDataJBAsp_fw2 <- AllDataJBAsp_fw2[!c(AllDataJBAsp_fw2$PFAA == "PFUA" ),]

# DataList$PFAA<-factor(DataList$PFAA)


## ----jba spring model bias tables, echo = FALSE, message=FALSE, warning=FALSE---------------------------

jba_sp1_error <- estimate_error(MetricData = AllDataJBAsp, DataID = "jba_sp1")
jba_sp1_error_fw <- estimate_error(MetricData = AllDataJBAsp_fw, DataID = "jba_sp1")

jba_sp1_error2 <- estimate_error(MetricData = AllDataJBAsp2, DataID = "jba_sp1")
jba_sp1_error_fw2 <- estimate_error(MetricData = AllDataJBAsp_fw2, DataID = "jba_sp1")



## ----jba spring figures, fig.cap= "Figure. 1", fig.height = 10, fig.width = 10.5------------------------

p1<-ggplot(data = AllDataJBAsp, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo, ymax = Obs_logngkg + Obs_logngkg_Up), size = 0.5) +
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo, xmax = log_ngkg_re + logngkg_re_Up), size = 0.5) +
    geom_point(pch = 19, size = 3)+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Known PFAS in diet')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Measured log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(2,6.5), xlim = c(2,6.5))+
    theme(legend.position = "none")
# p1_fish

p2<-ggplot(data = AllDataJBAsp_fw, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo, ymax = Obs_logngkg + Obs_logngkg_Up), size = 0.5) +
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo, xmax = log_ngkg_re + logngkg_re_Up), size = 0.5) +
    geom_point(pch = 19, size = 3)+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Modeled trophic transfer of PFAA')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Measured log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(0,6.5), xlim = c(0,6.5))+
    # theme(legend.position = "right")
    theme(legend.position = "none")
# p1_fish_re


# legend<-cowplot::get_legend(legend_plot)
# fig1<-cowplot::plot_grid(p1, p2, legend,
#                    ncol = 3,
#                    nrow = 1,
#                    rel_widths = c(1,1, 0.3),
#         labels = c("A", "B", ""))
# fig1

p1.2<-ggplot(data = AllDataJBAsp, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo, ymax = Obs_logngkg + Obs_logngkg_Up), size = 0.5) +
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo, xmax = log_ngkg_re + logngkg_re_Up), size = 0.5) +
    geom_point(pch = 19, size = 3)+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Known PFAS in diet: SAME DAY')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Measured log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(2,6.5), xlim = c(2,6.5))+
    theme(legend.position = "none")
# p1_fish

p2.2<-ggplot(data = AllDataJBAsp_fw, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo, ymax = Obs_logngkg + Obs_logngkg_Up), size = 0.5) +
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo, xmax = log_ngkg_re + logngkg_re_Up), size = 0.5) +
    geom_point(pch = 19, size = 3)+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Modeled trophic transfer of PFAA: SAME DAY')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Measured log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(0,6.5), xlim = c(0,6.5))+
    # theme(legend.position = "right")
    theme(legend.position = "none")
# p1_fish_re
legend_plot <- ggplot(data = AllDataJBAsp_fw,
         aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
  geom_point()+
  scale_color_manual(values = PFAA_cols)+
  theme(legend.title = element_blank())+
  theme_bw()

legend<-cowplot::get_legend(legend_plot)
fig1<-cowplot::plot_grid(p1, p2, legend,
                         p1.2, p2.2, NULL,
                   ncol = 3,
                   nrow = 2,
                   rel_widths = c(1,1, 0.3),
        labels = c("A", "B", "", "D", "E", ""))
fig1



## ----jba summer all data, echo = FALSE, message=FALSE, warning=FALSE------------------------------------
source(here("Code", "R data tables", "JBA_summer.R"))

settings_obsDiet = new_Settings(
  chooseDiet = 'forced'
)

# used below for sediment water ratio estimates
JBAsum_inputFiles_list<-inputFiles_list
water_JBAsum<-d.water %>%  # get the medians from the data
  select(PFHxS.m, PFOS.m, PFOA.m, PFNA.m) %>% 
  as.data.frame()
water_JBAsum<-water_JBAsum/1000 # original data in ng/L,  conver tto ng/mL
colnames(water_JBAsum)<-c("PFHxS", "PFOS", "PFOA", "PFNA") 
sed_JBAsum<-d.sed %>%  # get the medians from the data
  select(PFHxS.m, PFOS.m, PFOA.m, PFNA.m) %>% 
  as.data.frame()
colnames(sed_JBAsum)<-c("PFHxS", "PFOS", "PFOA", "PFNA") 
obs_ratio_JBAsum<-sed_JBAsum/water_JBAsum

AllDataJBAsum<-runEcosystemModel(settings = settings_obsDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR',
                                    "Gill_uptake","Dietary_uptake", "Sediment_uptake",
                                    "Gill_up%", "Diet_up%", 'Sed_up%',
                                    "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = inputFiles_list$max_dietData)


settings_modelDiet = new_Settings(
  chooseDiet = 'default'
)

AllDataJBAsum_fw<-runEcosystemModel(settings = settings_modelDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR',
                                    "Gill_uptake","Dietary_uptake", "Sediment_uptake",
                                    "Gill_up%", "Diet_up%", 'Sed_up%',
                                    "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = inputFiles_list$max_dietData)

AllDataJBAsum <- AllDataJBAsum[!c(AllDataJBAsum$PFAA == "PFUA" ),]
AllDataJBAsum_fw <- AllDataJBAsum_fw[!c(AllDataJBAsum_fw$PFAA == "PFUA" ),]


## ----jba summer same day data, echo = FALSE, message=FALSE, warning=FALSE-------------------------------
source(here("Code", "R data tables", "JBA_summer_sameDate.R"))

settings_obsDiet = new_Settings(
  chooseDiet = 'forced'
)

AllDataJBAsum2<-runEcosystemModel(settings = settings_obsDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR',
                                    "Gill_uptake","Dietary_uptake", "Sediment_uptake",
                                    "Gill_up%", "Diet_up%", 'Sed_up%',
                                    "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = inputFiles_list$max_dietData)


settings_modelDiet = new_Settings(
  chooseDiet = 'default'
)

AllDataJBAsum_fw2<-runEcosystemModel(settings = settings_modelDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR',
                                    "Gill_uptake","Dietary_uptake", "Sediment_uptake",
                                    "Gill_up%", "Diet_up%", 'Sed_up%',
                                    "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = inputFiles_list$max_dietData)

AllDataJBAsum2 <- AllDataJBAsum2[!c(AllDataJBAsum2$PFAA == "PFUA" ),]
AllDataJBAsum_fw2 <- AllDataJBAsum_fw2[!c(AllDataJBAsum_fw2$PFAA == "PFUA" ),]


## ----jba summer model bias tables, echo = FALSE, message=FALSE, warning=FALSE---------------------------

jba_sum1_error <- estimate_error(MetricData = AllDataJBAsum, DataID = "jba_sum1")
jba_sum1_error_fw <- estimate_error(MetricData = AllDataJBAsum_fw, DataID = "jba_sum1")

jba_sum1_error2 <- estimate_error(MetricData = AllDataJBAsum2, DataID = "jba_sum1")
jba_sum1_error_fw2 <- estimate_error(MetricData = AllDataJBAsum_fw2, DataID = "jba_sum1")



## ----jba plots and outputs, echo = FALSE, message=FALSE, warning=FALSE, fig.cap= "Figure 3", fig.height = 10, fig.width = 10.5----

p1<-ggplot(data = AllDataJBAsum, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo, ymax = Obs_logngkg + Obs_logngkg_Up), size = 0.5) +
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo, xmax = log_ngkg_re + logngkg_re_Up), size = 0.5) +
    geom_point(pch = 19, size = 3)+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Known PFAS in diet')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Measured log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(2,6.5), xlim = c(2,6.5))+
    theme(legend.position = "none")
# p1_fish

p2<-ggplot(data = AllDataJBAsum_fw, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo, ymax = Obs_logngkg + Obs_logngkg_Up), size = 0.5) +
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo, xmax = log_ngkg_re + logngkg_re_Up), size = 0.5) +
    geom_point(pch = 19, size = 3)+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Modeled trophic transfer of PFAA')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Measured log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(0,6.5), xlim = c(0,6.5))+
    # theme(legend.position = "right")
    theme(legend.position = "none")

p1.2<-ggplot(data = AllDataJBAsum, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo, ymax = Obs_logngkg + Obs_logngkg_Up), size = 0.5) +
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo, xmax = log_ngkg_re + logngkg_re_Up), size = 0.5) +
    geom_point(pch = 19, size = 3)+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Known PFAS in diet: same day')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Measured log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(2,6.5), xlim = c(2,6.5))+
    theme(legend.position = "none")
# p1_fish

p2.2<-ggplot(data = AllDataJBAsum_fw, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo, ymax = Obs_logngkg + Obs_logngkg_Up), size = 0.5) +
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo, xmax = log_ngkg_re + logngkg_re_Up), size = 0.5) +
    geom_point(pch = 19, size = 3)+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Modeled trophic transfer of PFAA: same day')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Measured log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(0,6.5), xlim = c(0,6.5))+
    # theme(legend.position = "right")
    theme(legend.position = "none")


# p1_fish_re
legend_plot <- ggplot(data = AllDataJBAsum_fw,
         aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
  geom_point()+
  scale_color_manual(values = PFAA_cols)+
  theme(legend.title = element_blank())+
  theme_bw()

legend<-cowplot::get_legend(legend_plot)
# fig1<-cowplot::plot_grid(p1, p2, legend,
#                    ncol = 3,
#                    nrow = 1,
#                    rel_widths = c(1,1, 0.3),
#         labels = c("A", "B", ""))
# fig1

fig1<-cowplot::plot_grid(p1, p2, legend,
                         p1.2, p2.2, NULL,
                   ncol = 3,
                   nrow = 2,
                   rel_widths = c(1,1, 0.3),
        labels = c("A", "B", "", "D", "E", ""))
fig1



## ----setup2---------------------------------------------------------------------------------------------
AllDataJBAsp$System<-"JBA-spring"
AllDataJBAsp_fw$System<-"JBA-spring"

AllDataJBAsum$System<-"JBA-summer"
AllDataJBAsum_fw$System<-"JBA-summer"

AllDataWG$System<-"Willow Grove"
AllDataWG_fw$System<-"Willow Grove"

d<-rbind(AllDataJBAsp, AllDataJBAsum, AllDataWG)
d_fw<-rbind(AllDataJBAsp_fw, AllDataJBAsum_fw, AllDataWG_fw)
d$System<-as.factor(d$System)
d_fw$System<-as.factor(d_fw$System)

p1<-ggplot(data = d, aes(x = log_ngkg_re, y = Obs_logngkg,
                         color = PFAA,
                         shape = System))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y, shape = NULL),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1, shape = NULL),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2, shape = NULL),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3, shape = NULL),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4, shape = NULL),
              linetype = 3, linewidth = 0.2, color = "black")+
    # geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo,
    #                   ymax = Obs_logngkg + Obs_logngkg_Up), size = 0.5) +
    # geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo,
    #                    xmax = log_ngkg_re + logngkg_re_Up), size = 0.5) +
    geom_text(mapping = aes(label = SppAlias), check_overlap = T, size = 2)+
    # geom_point(size = 3)+
    facet_grid(.~System)+
    scale_shape_manual(values = c(19, 21, 23))+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Given PFAA diet values')+
    coord_cartesian(ylim = c(1,7), xlim = c(1,7))

ggformat(p1, y_title = ('Measured log PFAS (ng/kg)'), x_title = 'Modeled log PFAS (ng/kg)', size_text = 12, print = F)
p1<-p1+theme(legend.position = "none")

p2<-ggplot(data = d_fw, aes(x = log_ngkg_re, y = Obs_logngkg,
                            color = PFAA,
                            shape = System))+
    geom_line(data = p_lines, mapping = aes(x = x, y = y, shape = NULL),
              linetype = 1, linewidth = 0.5, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y1, shape = NULL),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y2, shape = NULL),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y3, shape = NULL),
              linetype = 3, linewidth = 0.2, color = "black")+
    geom_line(data = p_lines, mapping = aes(x = x, y = y4, shape = NULL),
              linetype = 3, linewidth = 0.2, color = "black")+
    # geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo,
    #                   ymax = Obs_logngkg + Obs_logngkg_Up), size = 0.5) +
    # geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_re_Lo,
    #                    xmax = log_ngkg_re + logngkg_re_Up), size = 0.5) +
    geom_text(mapping = aes(label = SppAlias), check_overlap = T, size = 2)+
    # geom_point(size = 3)+
    facet_grid(.~System)+
    scale_shape_manual(values = c(19, 21, 23))+
    theme_bw()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Observed diet')+
    coord_cartesian(ylim = c(1,7), xlim = c(1,7))

ggformat(p2, y_title = ('Measured log PFAS (ng/kg)'), x_title = 'Modeled log PFAS (ng/kg)', size_text = 12, print = F)
p2<-p2+theme(legend.position = "none")

# p1_fish_re
legend_plot <- ggplot(data = d,
         aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA, shape = System))+
  geom_point()+
  scale_color_manual(values = PFAA_cols)+
  scale_shape_manual(values = c(19, 21, 23))+
  theme(legend.title = element_blank())+
  theme_bw()


## ----legend and saving, fig.width=10, fig.height=7, fig.cap= "Figure 4"---------------------------------
legend<-cowplot::get_legend(legend_plot)
fig1<-cowplot::plot_grid(p1, p2, legend,
                   ncol = 3,
                   nrow = 1,
                   rel_widths = c(0.9, 0.9, 0.3),
                   label_x = c(0, 0, 0),
                   label_y = c(1,1,1),
                   label_size = 12,
        labels = c("known PFAS in diet", "modeled PFAS in diet", "")
        )
fig1

# ggsave(filename = "../../Figures/fig1_oct31.png", plot = fig1, height = 3.5, width = 8)



## ----error estimates MB all systems together------------------------------------------------------------
jba_sum1_error_pfaa<-as.data.frame(jba_sum1_error[[1]][,]) 
jba_sum1_error_pfaa$System<-"JBA-summer"
jba_sum1_error_pfaa_fw<-as.data.frame(jba_sum1_error_fw[[1]][,])
jba_sum1_error_pfaa_fw$System <- "JBA-summer"
jba_sp1_error_pfaa<-as.data.frame(jba_sp1_error[[1]][,])
jba_sp1_error_pfaa$System <- "JBA-spring"
jba_sp1_error_pfaa_fw<-as.data.frame(jba_sp1_error_fw[[1]][,])
jba_sp1_error_pfaa_fw$System <- "JBA-spring"
wg1_error_pfaa<-as.data.frame(wg1_error[[1]][,])
wg1_error_pfaa$System <- "Willow Grove"
wg1_error_pfaa_fw<-as.data.frame(wg1_error_fw[[1]][,])
wg1_error_pfaa_fw$System <- "Willow Grove"

colnames(jba_sum1_error_pfaa)<-c("PFAA", "MB", "System")
colnames(jba_sum1_error_pfaa_fw)<-c("PFAA", "MB", "System")
colnames(jba_sp1_error_pfaa)<-c("PFAA", "MB", "System")
colnames(jba_sp1_error_pfaa_fw)<-c("PFAA", "MB", "System")
colnames(wg1_error_pfaa)<-c("PFAA", "MB", "System")
colnames(wg1_error_pfaa_fw)<-c("PFAA", "MB", "System")

jba_sum1_error_sp<-as.data.frame(jba_sum1_error[[2]][,])
jba_sum1_error_sp$System <- "JBA-summer"
jba_sum1_error_sp_fw<-as.data.frame(jba_sum1_error_fw[[2]][,])
jba_sum1_error_sp_fw$System <- "JBA-summer"
jba_sp1_error_sp<-as.data.frame(jba_sp1_error[[2]][,])
jba_sp1_error_sp$System <- "JBA-spring"
jba_sp1_error_sp_fw<-as.data.frame(jba_sp1_error_fw[[2]][,])
jba_sp1_error_sp_fw$System <- "JBA-spring"
wg1_error_sp<-as.data.frame(wg1_error[[2]][,])
wg1_error_sp$System <- "Willow Grove"
wg1_error_sp_fw<-as.data.frame(wg1_error_fw[[2]][,])
wg1_error_sp_fw$System <- "Willow Grove"

colnames(jba_sum1_error_sp)<-c("Species", "MB", "System")
colnames(jba_sum1_error_sp_fw)<-c("Species", "MB", "System")
colnames(jba_sp1_error_sp)<-c("Species", "MB", "System")
colnames(jba_sp1_error_sp_fw)<-c("Species", "MB", "System")
colnames(wg1_error_sp)<-c("Species", "MB", "System")
colnames(wg1_error_sp_fw)<-c("Species", "MB", "System")

d_err_pfaa<-rbind(jba_sum1_error_pfaa, jba_sp1_error_pfaa, wg1_error_pfaa)
d_err_pfaa_fw<-rbind(jba_sum1_error_pfaa_fw, jba_sp1_error_pfaa_fw, wg1_error_pfaa_fw)
d_err_pfaa$model<-"known diet"
d_err_pfaa_fw$model<-"modeled diet"
d_err_sp<-rbind(jba_sum1_error_sp, jba_sp1_error_sp, wg1_error_sp)
d_err_sp_fw<-rbind(jba_sum1_error_sp_fw, jba_sp1_error_sp_fw, wg1_error_sp_fw)
d_err_sp$model<-"known diet"
d_err_sp_fw$model<-"modeled diet"



## ----model bias plots, fig.cap= "Figure 5."-------------------------------------------------------------
p_mb1<-ggplot(data = rbind(d_err_pfaa, d_err_pfaa_fw),
           aes(color = PFAA, 
               x = System, y = MB,
               shape = model, 
               group = interaction(model, PFAA)))+
  geom_point(position = position_dodge(width = 0.5),
           stat = "identity", size=3)+
  scale_shape_manual(values = c(21,19))+
  geom_hline(yintercept = 1, color = "black")+
  geom_hline(yintercept = c(0.5, 2), color = "black", linetype = 2, linewidth = 0.1)+
  # geom_text(mapping = aes(group = PFAA, vjust = 0.5 - sign(value)/1), size = 3,
  #           color = "black", position = position_dodge(width = 0.8))+
  scale_color_manual(values = PFAA_cols)+
  # scale_fill_manual(values = PFAA_cols)+
  theme_bw()+
  # ylim(-15,15)+
  theme(legend.position = "top")+
  ylab("Model Bias (tissue conc.)")+
  theme(axis.text.x = element_text(angle = 0))+
  coord_flip()+
  scale_x_discrete(labels = c("JBA \n spring", "JBA \n summer", "Willow \nGrove"))

# ggformat(p_mb1, y_title = "Model Bias (tissue conc.)", x_title = "", size_text = 12, print = T)
# p_mb1<- p_mb1 +
#   theme(
#     panel.background = element_rect(fill='transparent'), #transparent panel bg
#     plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
#     panel.grid.major = element_blank(), #remove major gridlines
#     panel.grid.minor = element_blank(), #remove minor gridlines
#     legend.background = element_rect(fill='transparent'), #transparent legend bg
#     legend.box.background = element_rect(fill='transparent') #transparent legend panel
#   )
p_mb1

# 
# ggsave(filename = "../../Figures/fig2_mb_oct31.png", plot = p_mb1, height = 3.5, width = 6, bg = "transparent")




## ----report in poster stats-----------------------------------------------------------------------------

pfaa_err<-rbind(d_err_pfaa, d_err_pfaa_fw)
pfaa_err$AMB<-abs(pfaa_err$MB)
n_2<-nrow(pfaa_err[pfaa_err$AMB <= 2, ])
n_10<-nrow(pfaa_err[pfaa_err$AMB <= 10, ])
n_2/nrow(pfaa_err)
n_10/nrow(pfaa_err)



## ----model bias plots 2, fig.cap= "Figure 6."-----------------------------------------------------------

p_mb2<-ggplot(data = rbind(d_err_sp, d_err_sp_fw),
           aes(color = System, 
               x = Species, y = MB,
               shape = model, 
               group = interaction( System, model)))+
  geom_segment(aes(x=Species, xend = Species, y=1, yend=MB, 
                   color = System,
               group = interaction( System, model)),linewidth = 3, alpha= 0.7) +
  scale_shape_manual(values = c(19,19))+
  geom_point(stat = "identity", size = 3)+
  geom_hline(yintercept = 1, color = "black")+
  geom_hline(yintercept = c(2, 0.5), color = "grey30", linetype = 2, linewidth = 0.1)+
  # geom_text(mapping = aes(group = PFAA, vjust = 0.5 - sign(value)/1), size = 3,
  #           color = "black", position = position_dodge(width = 0.8))+
  scale_color_manual(values = c("black", "grey50", "green4"))+
  scale_fill_manual(values = c("black", "grey50", "green4"))+
  theme_bw()+
  theme(legend.position = "none")+
  ylab("Model Bias (tissue conc.)")+
  theme(axis.text.x = element_text(angle = 0))+
  coord_flip()+
  facet_grid(model~System)
  # scale_x_discrete(labels = c("JBA \n spring", "JBA \n summer", "Willow \nGrove"))
p_mb2

# ggformat(p_mb2, y_title = "Model Bias (tissue conc.)", x_title = "", size_text = 12, print = T)
# ggsave(filename = "../../Figures/fig2_mb2_oct31.png", plot = p_mb2, height = 3.5, width = 6)


p_mb1<-ggplot(data = rbind(d_err_pfaa, d_err_pfaa_fw),
           aes(color = System, 
               x = PFAA, y = MB,
               shape = model, 
               group = interaction( System, model)))+
  geom_segment(aes(x=PFAA, xend = PFAA, y=1, yend=MB, 
                   color = System,
               group = interaction( System, model)),linewidth = 3, alpha= 0.7) +
  scale_shape_manual(values = c(19,19))+
  geom_point(stat = "identity", size = 3)+
  geom_hline(yintercept = 1, color = "black")+
  geom_hline(yintercept = c(2, 0.5), color = "grey30", linetype = 2, linewidth = 0.1)+
  # geom_text(mapping = aes(group = PFAA, vjust = 0.5 - sign(value)/1), size = 3,
  #           color = "black", position = position_dodge(width = 0.8))+
  scale_color_manual(values = c("black", "grey50", "green4"))+
  scale_fill_manual(values = c("black", "grey50", "green4"))+
  theme_bw()+
  theme(legend.position = "none")+
  ylab("Model Bias (tissue conc.)")+
  theme(axis.text.x = element_text(angle = 0))+
  coord_flip()+
  facet_grid(model~System)
  # scale_x_discrete(labels = c("JBA \n spring", "JBA \n summer", "Willow \nGrove"))
p_mb1



## ----setup dataframes add trophic levels----------------------------------------------------------------

# JBA summer
jba_sum_trophic_levels<-data.frame(
  SppAlias = c("Phy", "Dac", "Dar", "Min", "Fal", "Bas", "Mad", "Pum", "Swa"), 
  trophic = c(0,       2.62, 3.75,  2.92,  3.25,  1.5,   3.53,  3.83,  2.5),
  System = c("JBA-summer"),
  loc = "JBA sum")

jba_sp_trophic_levels<-data.frame(
  SppAlias = c("Phy", "Kil", "Chu", "Dac", "Dar", "Min", "Mad", "Pum", "Swa"), 
  trophic = c(0,       2.5,  3.93,   2.62,  3.75,  2.92,  3.53,  3.83,  2.5),
  System = c("JBA-spring"),
  loc = "JBA sp")

wg_trophic_levels<-data.frame(
  SppAlias = c("Phy", "Pry", "Bgl", "Bas"), 
  trophic = c(0,        2.31, 3.08,  3.69),
  System = c("Willow Grove"), 
  loc = "WG")

d_sed<-merge(d,
      rbind(jba_sum_trophic_levels, jba_sp_trophic_levels, wg_trophic_levels),
      by = c("SppAlias", "System"))

AllDataWG$loc<-"WG"
AllDataJBAsp$loc<-"JBA sp"
AllDataJBAsum$loc<-"JBA sum"
AllData<-rbind(AllDataWG, AllDataJBAsp, AllDataJBAsum) # only 
AllData <- AllData[!c(AllData$PFAA == "PFUA" ),]

AllData<-merge(AllData,
      rbind(jba_sum_trophic_levels, jba_sp_trophic_levels, wg_trophic_levels),
      by = c("SppAlias", "loc", "System"))
      

# stackbar dataset
AllData_m<-melt(AllData[, c("SppAlias", "PFAA", "Gill_up%", "Diet_up%", "Sed_up%",
                            "WB", "C_WTO", "C_s", "loc")],
                id.vars=c("loc","SppAlias", "PFAA", "WB", "C_WTO", "C_s"))

# summarize means for % uptake for each run
mean_d_pct<-AllData_m %>%
dplyr::group_by(PFAA, SppAlias, variable, C_WTO, C_s, WB, loc) %>%
summarize(mean_pct = mean(value),
          sd_pct = sd(value)) 

mean_d_pct_m<-mean_d_pct %>% 
  dplyr::group_by(PFAA, variable, C_WTO, C_s) %>%
  summarize(mean_pct = mean(mean_pct),
          sd_pct = sd(mean_pct)) 
  
# summarize actual values 
mean_d<-AllData %>%
dplyr::group_by(PFAA, SppAlias, loc) %>%
summarize(Gill = mean(`Gill_uptake`),
          Diet = mean(`Dietary_uptake`),
          Sediment = mean(`Sediment_uptake`),
          Elimination = mean(`Total Elimination`), 
          WB = mean(WB)) %>%
pivot_longer(cols = c(Gill, Diet, Sediment, Elimination), names_to = "Path", values_to = "Uptake_g.kg.day")

mean_sd<-AllData %>%
dplyr::group_by(PFAA, SppAlias, loc) %>%
summarize(Gill = sd(`Gill_uptake`),
          Diet = sd(`Dietary_uptake`),
          Sediment = sd(`Sediment_uptake`),
          Elimination = sd(`Total Elimination`),
          WB = mean(WB)) %>%
pivot_longer(cols = c(Gill, Diet, Sediment, Elimination), names_to = "Path", values_to = "Uptake_SD_g.kg.day")

AllData_m2 <- merge(mean_d, mean_sd, all.x = TRUE, by = c("PFAA", "SppAlias", "Path", "loc", "WB"))

# water to sediment PFAS ratio
AllData$sed_wt_pfas_ratio <- AllData$C_s/AllData$C_WTO

# labels chemID
PFAA_List <- data.frame(PFAA = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'), PFAA_labels = c('C6 PFSA','C8 PFSA','C8 PFCA', 'C9 PFCA','C10 PFCA','C11 PFCA'), PFAA_n = c(6,8,8,9,10,11))
AllData<-merge(AllData, PFAA_List, by = "PFAA")

AllData$PFAA_labels<-factor(AllData$PFAA_labels)
AllData$PFAA_labels<-factor(AllData$PFAA_labels, levels = c('C6 PFSA','C8 PFSA','C8 PFCA', 'C9 PFCA','C10 PFCA'))


## ----plots sediment and error, eval = FALSE-------------------------------------------------------------
## 
## ggplot(d_sed[!d_sed$SppAlias=="Phy",],
##        aes(y = Dietary_uptake,
##            x = trophic,
##               color = PFAA,
##               shape = System))+
##   geom_point()+
##   facet_grid(.~System, scales = "free")+
##   # geom_point( mapping = aes(y = Obs_logngkg, x = Sed_trophic,
##   #               color = PFAA,
##   #               shape = System),
##   #            position = position_dodge(width = 0.3, preserve = "single"),
##   #            size = 2, alpha = 0.4)+
##   scale_color_manual(values = PFAA_cols)+
##   scale_shape_manual(values = c(19, 21, 23))+
##   theme_bw()+
##   geom_smooth(se = FALSE, method = "lm", formula=y ~ poly(x, 1, raw=TRUE) )+
##   xlab("Trophic level")+
##   xlab("Dietary uptake of PFAS (ng/kg/d)")+
##   theme(axis.text.x = element_text(angle = 0))
## 
## ggformat(p_mb1, y_title = "Model Bias (tissue conc.)",
##          x_title = "",
##          size_text = 12, print = T)
## 
## # ggsave(filename = "../../Figures/fig2_mb_oct31.png", plot = p_mb1, height = 3.5, width = 6)
## 
## 


## ----water vs gill figure, echo = FALSE, message = FALSE, warning  = FALSE, fig.height = 10, fig.width = 10, fig.align = "center", eval = TRUE, fig.cap= "Figure 7."----
ggplot(mean_d_pct[!mean_d_pct$SppAlias == "Phy", ],
       aes( x=mean_pct*100, y=SppAlias, color = variable, fill = variable)) + 
    geom_bar(stat="identity", position = "stack") +
    # geom_errorbarh(data = mean_d_pct[mean_d_pct$variable == "Diet_up%" & !mean_d_pct$SppAlias == "Phy", ],
    #                 aes(xmin=mean_pct*100-sd_pct*100, xmax=mean_pct*100+sd_pct*100), 
    #                 position = "dodge", size = 0.5, height = 0.3, color = "black")+
    facet_wrap(loc~PFAA, scales = "free")+
    scale_fill_manual(breaks = c("Diet_up%", "Gill_up%", "Sed_up%"), values = c("#007445", "#00C6FF", "#915801"),
                      labels=c("Diet uptake","Gill uptake", "Sediment uptake"))+
    scale_color_manual(breaks = c("Diet_up%", "Gill_up%", "Sed_up%"), values = c("#007445", "#00C2FF", "#915801"),
                      labels=c("Diet uptake","Gill uptake", "Sediment uptake"))+  ylab("Species")+
  xlab("% Uptake")+
  # coord_polar()+
  theme_bw()+
  theme(legend.position = "top")


## ----diet to water uptake, fig.width=5, fig.height=10, eval = FALSE-------------------------------------
## 
## p_diet_water<-ggplot(mean_d_pct[!mean_d_pct$SppAlias == "Phy" & mean_d_pct$loc == "JBA sum", ],
##        aes( x=mean_pct*100, y=SppAlias,
##             color = variable,
##             fill = variable,
##             group = loc)) +
##     geom_bar(stat="identity", position = "stack",
##              color = "black",
##              stroke = 0.3)+
##     facet_wrap(.~PFAA, scales = "free", nrow =5)+
##     scale_fill_manual(breaks = c("Diet_up%", "Gill_up%"), values = c("#007445", "#00C2FF", "#915801"),
##                       labels=c("Diet uptake","Gill uptake"))+
##   scale_y_discrete(breaks = c("Swa",
##                               "Pum",
##                               "Min",
##                               "Mad",
##                               "Fal",
##                               "Dar",
##                               "Dac",
##                               "Bas"),
##                    labels = c("Swallowtail Shiner",
##                               "Pumkinseed",
##                               "Mudminnow",
##                               "Margined Madtom",
##                               "Fallfish",
##                               "Darter sp.",
##                               "Dace sp.",
##                               "Largemouth Bass"))+
##   ylab("Species")+
##   xlab("% Uptake")+
##   theme_bw()+
##   theme(legend.position = "right")
## 
## ggformat(p_diet_water, y_title = 'Species', x_title = "% Uptake", title = "JBA summer", size_text = 12)
## p_diet_water<-p_diet_water+ theme(legend.position = "right")
## # ggsave(filename = "../Figures/fig8_deit_water_now1.png", plot = p_diet_water, height = 10, width = 5)
## 
## 


## ----uptake and depuration rates per day, echo = FALSE, message = FALSE, warning  = FALSE, fig.width=10, fig.height=7,  fig.align = "center", eval = FALSE, fig.cap= "Figure 8."----
## ggplot(AllData_m2[!AllData_m2$SppAlias == "Phy", ],
##       aes(fill = Path, x = Uptake_g.kg.day, y = SppAlias, color = Path)) +
##     geom_bar(stat="identity", position = "dodge", color = "white") +
##     # geom_pointrange(aes(y = SppAlias,
##     #                xmin=`Uptake_%g.kg.day`-`Uptake_SD_%g.kg.day`,
##     #                xmax=`Uptake_%g.kg.day`+`Uptake_SD_%g.kg.day`), position = position_dodge(width = 1))+
##     facet_wrap(.~PFAA, scales = "free")+
##     scale_fill_manual(breaks = c("Diet", "Gill", "Sediment", "Elimination"),
##                       values = c("#007445", "#00C2FF", "#915801", "black"),
##                       labels=c("Diet uptake","Gill uptake", "Sediment", "Total Elimination"))+
##   ylab("Species")+
##   xlab("Uptake h/kg/d")+
##   theme_bw()+
##   theme(legend.position = "top",
##         legend.title = element_blank(),
##         legend.text=element_text(size=15))
## 
## 


## ----Diet and total tissue PFAA-------------------------------------------------------------------------
# ggplot(data = AllData[c(!AllData$SppAlias =="Phy" & AllData$loc == "JBA sum"),],
#        aes(y = log_ngkg_re,
#            color = PFAA, x = sed_wt_pfas_ratio, 
#            group  = interaction(PFAA, SppAlias, loc)))+
#   geom_smooth(size = 0.05, se = F)+
#   geom_point(pch=19)+
#   scale_color_manual(values = PFAA_cols)


pdiet_water2<-ggplot(data = AllData[c(!AllData$SppAlias =="Phy" ),],
       aes(x = `Diet_up%`, 
           y = log_ngkg_re, 
           color = PFAA,
           shape = loc,
           group  = interaction(PFAA)))+
  geom_point(size = 3)+
  scale_color_manual(values = PFAA_cols)

ggformat(pdiet_water2, x_title = ("Dietary uptake (% total)"), y_title = ("Modeled log PFAS (ng/kg); known diet"), size_text = 12, print = T)




## ----data set wrangle organise--------------------------------------------------------------------------
# stackbar dataset
AllData_m<-melt(AllData[, c("SppAlias", "PFAA", "Gill_up%", "Diet_up%", "WB", "C_WTO", "C_s", "loc")], id.vars=c("loc","SppAlias", "PFAA", "WB", "C_WTO", "C_s"))

# summarize means for % uptake for each run
mean_d_pct<-AllData_m %>%
dplyr::group_by(PFAA, SppAlias, variable, C_WTO, C_s, WB, loc) %>%
summarize(mean_pct = mean(value),
          sd_pct = sd(value)) 

mean_d_pct_m<-mean_d_pct %>% 
  dplyr::group_by(PFAA, variable, C_WTO, C_s) %>%
  summarize(mean_pct = mean(mean_pct),
          sd_pct = sd(mean_pct)) 
  
# summarize actual values 
mean_d<-AllData %>%
dplyr::group_by(PFAA, SppAlias, loc) %>%
summarize(Gill = mean(`Gill_uptake`),
          Diet = mean(`Dietary_uptake`),
          Elimination = mean(`Total Elimination`), 
          WB = mean(WB)) %>%
pivot_longer(cols = c(Gill, Diet, Elimination), names_to = "Path", values_to = "Uptake_g.kg.day")

mean_sd<-AllData %>%
dplyr::group_by(PFAA, SppAlias, loc) %>%
summarize(Gill = sd(`Gill_uptake`),
          Diet = sd(`Dietary_uptake`),
          Elimination = sd(`Total Elimination`),
          WB = mean(WB)) %>%
pivot_longer(cols = c(Gill, Diet, Elimination), names_to = "Path", values_to = "Uptake_SD_g.kg.day")

AllData_m2 <- merge(mean_d, mean_sd, all.x = TRUE, by = c("PFAA", "SppAlias", "Path", "loc", "WB"))

# water to sediment PFAS ratio
AllData$sed_wt_pfas_ratio <- AllData$C_s/AllData$C_WTO



## ----trophic level and the effect on the diet, fig.height = 6, fig.width = 8, fig.align = "center", eval = TRUE----
ggplot(AllData[!AllData$SppAlias == "Phy",], aes(y=`Diet_up%`*100,
                    x = trophic, color = PFAA))+
  geom_point(size = 2)+
  facet_grid(.~loc, scales = "free")+
  scale_color_manual(values = PFAA_cols)+
  # ylim(0,100)+
  geom_smooth(method = "lm", se = F, linewidth = 1.5, alpha = 0.5)+
  theme_bw()+
  xlab("Trophic level")+
  ylab("PFAA dietary uptake (%)")
# ggsave(filename = "./Figures/poster_sediment_PFOS_bioaccum_JBAsummer.png", width = 5, height = 5)


## ----Sed_trophic and error rates------------------------------------------------------------------------
ggplot(AllData[!AllData$SppAlias == "Phy",],
       aes(y=`Diet_up%`*100,
           x = PFAA_labels,
           group = SppAlias,
           color = PFAA))+
  geom_point(size = 1)+
  facet_wrap(.~loc)+
  scale_color_manual(values = PFAA_cols)+
  # geom_line(linewidth = 0.5)+
  ylim(0,100)+
  theme_bw()+
  xlab("PFAA chain length")+
  ylab("PFAA dietary uptake (%)")+
  theme(axis.text.x = element_text(angle = 45))



## ----Sed_trophic and error rates2-----------------------------------------------------------------------
ggplot(AllData,
       aes(y= log_ngkg_re,
           x = PFAA_labels,
           group = SppAlias, 
           color = PFAA))+
  geom_point(size = 1)+
  geom_point(size = 1, pch=21, position = position_nudge(x = 0.1),
             mapping = aes(y= Obs_logngkg,
                           x = PFAA_labels,
                           group = SppAlias, 
                           color = PFAA))+
  scale_color_manual(values = PFAA_cols)+
  # geom_line(linewidth = 0.5)+
  theme_bw()+
  xlab("PFAA chain length")+
  ylab("Tissue PFAA bioconcentriotions log(ng/kg)")+
  theme(axis.text.x = element_text(angle = 45))



## ----sediment to water ratio----------------------------------------------------------------------------
p_sed1<-ggplot(AllData[!c(AllData$SppAlias == "Phy"),],
       aes(y=log_ngkg_re, x = sed_wt_pfas_ratio,
           group = interaction( PFAA),
           color = PFAA, size = C_s))+
  # geom_smooth(method = "lm", se = F)+
  geom_point( pch=21)+
  # facet_grid(.~PFAA, scales = "free")+
  scale_color_manual(values = PFAA_cols)+
  scale_linetype_manual(values = c("dashed", "dotted", "solid"))+
  theme_bw()
ggformat(p_sed1, x_title = ("PFAA sediment:water ratio"),
         y_title = ("Tissue PFAA (ng/kg)"))


## ----sediment to water ratio2---------------------------------------------------------------------------

means_d<-AllData[!c(AllData$SppAlias == "Phy"),] %>% 
  dplyr:::group_by(loc, PFAA) %>% 
  summarize(mean_pfos = mean(`Diet_up%`), 
            mean_ratio = mean(sed_wt_pfas_ratio), 
            C_s = mean(C_s))
  

p_sed1<-ggplot(AllData[!c(AllData$SppAlias == "Phy"),],
       aes(y= `Diet_up%`, x = sed_wt_pfas_ratio,
           group = interaction( PFAA),
           color = PFAA,
           size = C_s))+
  # geom_point(data = means_d, aes(x = mean_ratio, y = mean_pfos, 
                                 # size = NULL), size = 4)+
  # geom_line(data = means_d, aes(x = mean_ratio, y = mean_pfos, 
                                # size = NULL), linewidth = 0.5)+
  # geom_smooth(method = "lm", formula = y~poly(x,1), se = F)+
  geom_point( pch=21)+
  # facet_grid(.~PFAA, scales = "free")+
  scale_color_manual(values = PFAA_cols)+
  scale_linetype_manual(values = c("dashed", "dotted", "solid"))+
  theme_bw()
ggformat(p_sed1, x_title = ("PFAA sediment:water ratio"),
         y_title = ("% Diet uptake"), print = T)




## ----loop sediment ratios source data, message = FALSE, eval=T, warning  = FALSE,  fig.width=10, fig.height=7,  fig.align = "center"----

# the data are read in and saved above when sourcing tables for each system
# create an identical data frame that will be filled in with values 
# based on the selected sed:water ratio
sed_WG<-water_WG[-1,] 
sed_JBAsp<-water_JBAsp[-1,]       
sed_JBAsum<-water_JBAsum[-1,]  

settings_modelDiet = new_Settings(
  chooseDiet = 'default'
)

# sediment water ratios that are empirically observed
l<-2
# r_PFAAs <- c(seq(0.01,0.1, length.out = l/2), seq(1, 5, length.out = l/2))
r_PFAAs<-c(seq(0.001,5, length.out = 0.5*l), seq(7, 30, length.out = 0.5*l))
# custom length out, for 30 low range, 30 high range. 

for(i in 1:l){
  # loop through all ratios 
  
  sed_WG[1,"PFHxS"] <- (water_WG[1, "PFHxS"]) * r_PFAAs[i]  # water in ng/ml, sed in ng/g
  sed_WG[1,"PFOS"] <- (water_WG[1, "PFOS"]) * r_PFAAs[i] 
  sed_WG[1,"PFOA"] <- (water_WG[1, "PFOA"]) * r_PFAAs[i] 
  sed_WG[1,"PFNA"] <- (water_WG[1, "PFNA"]) * r_PFAAs[i] 

  sed_JBAsp[1,"PFHxS"]<-(water_JBAsp[1, "PFHxS"]) * r_PFAAs[i] 
  sed_JBAsp[1,"PFOS"]<-(water_JBAsp[1, "PFOS"]) * r_PFAAs[i]
  sed_JBAsp[1,"PFOA"]<-(water_JBAsp[1, "PFOA"]) * r_PFAAs[i]
  sed_JBAsp[1,"PFNA"]<-(water_JBAsp[1, "PFNA"]) * r_PFAAs[i] 

  sed_JBAsum[1,"PFHxS"]<-(water_JBAsum[1, "PFHxS"]) * r_PFAAs[i] 
  sed_JBAsum[1,"PFOS"]<-(water_JBAsum[1, "PFOS"]) * r_PFAAs[i]
  sed_JBAsum[1,"PFOA"]<-(water_JBAsum[1, "PFOA"]) * r_PFAAs[i] 
  sed_JBAsum[1,"PFNA"]<-(water_JBAsum[1, "PFNA"]) * r_PFAAs[i] 

  WG_inputFiles_list$chemicalData["C_s", 1:4]<- sed_WG[1, c("PFHxS", "PFOS", "PFOA", "PFNA")]
  WG_inputFiles_list$chemicalData["C_WTO",1:4]<-water_WG [1, c("PFHxS", "PFOS", "PFOA", "PFNA")]
  JBAsp_inputFiles_list$chemicalData["C_s", 1:4]<- sed_WG[1, c("PFHxS", "PFOS", "PFOA", "PFNA")]
  JBAsp_inputFiles_list$chemicalData["C_WTO", 1:4]<-water_WG [1, c("PFHxS", "PFOS", "PFOA", "PFNA")]
  JBAsum_inputFiles_list$chemicalData["C_s", 1:4]<- sed_WG[1, c("PFHxS", "PFOS", "PFOA", "PFNA")]
  JBAsum_inputFiles_list$chemicalData["C_WTO", 1:4]<-water_WG [1, c("PFHxS", "PFOS", "PFOA", "PFNA")]
  
  # WG 
  AllData_wg0<-runEcosystemModel(settings = settings_modelDiet,
                  inputFiles_list = WG_inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR', "Gill_uptake",
                                    "Dietary_uptake", "Gill_up%", "Diet_up%",
                                    "Total Elimination", "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA'),
                  dietData = WG_inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = WG_inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = WG_inputFiles_list$max_dietData)


  # JBA spring
  AllData_jba_sp0<-runEcosystemModel(settings = settings_modelDiet,
                  inputFiles_list = JBAsp_inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR', "Gill_uptake",
                                    "Dietary_uptake", "Gill_up%", "Diet_up%",
                                    "Total Elimination", "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA'),
                  dietData = JBAsp_inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = JBAsp_inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = JBAsp_inputFiles_list$max_dietData)

  # JBA summer 
  AllData_jba_sum0<-runEcosystemModel(settings = settings_modelDiet,
                  inputFiles_list = JBAsum_inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR', "Gill_uptake",
                                    "Dietary_uptake", "Gill_up%", "Diet_up%",
                                    "Total Elimination", "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA'),
                  dietData = JBAsum_inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = TRUE,
                  min_dietData = JBAsum_inputFiles_list$min_dietData,
                  max_diet_WTO_C_s = TRUE, 
                  max_dietData = JBAsum_inputFiles_list$max_dietData)

  AllData_wg0$loc<-"WG"
  AllData_jba_sp0$loc<-"JBA sp"
  AllData_jba_sum0$loc<-"JBA sum"
  AllData_wg0$ratio<- r_PFAAs[i]
  AllData_jba_sp0$ratio<- r_PFAAs[i]
  AllData_jba_sum0$ratio<- r_PFAAs[i]
  AllData0<-rbind(AllData_wg0, AllData_jba_sum0, AllData_jba_sp0)
  AllData0$sed_wt_pfas_ratio<-AllData0$C_s/AllData0$C_WTO
  all(AllData0$sed_wt_pfas_ratio == AllData_jba_sp0$ratio)
  
  if(i == 1 ){
    AllData<-AllData0
  }else{
    AllData<-rbind(AllData, AllData0)
  }
  
}

AllData_sed<-AllData


## ----sequence sediment plots ,message = FALSE, warning  = FALSE,  fig.height = 3, fig.width = 10, fig.align = "center"----
p_sed1<-ggplot(AllData_sed[!c(AllData_sed$SppAlias == "Phy"),],
       aes(y=`Diet_up%`*100, x = C_s,
           group = interaction(loc, SppAlias, PFAA),
           color = PFAA))+
  # geom_point(size = 0.5)+
  facet_grid(.~PFAA, scales = "free")+
  scale_color_manual(values = PFAA_cols)+
  scale_linetype_manual(values = c("dashed", "dotted", "solid"))+
  ylim(0,100)+
  geom_line(size= 0.5)+
  theme_bw()+
  xlab("PFAA in sediment (ng/g)")+
  ylab("PFAA dietary uptake (%)")
p_sed1


## ----sediment plots2 ,message = FALSE, warning  = FALSE,  fig.height = 5, fig.width = 5, fig.align = "center"----
p_sed2<-ggplot(AllData_sed[!c(AllData_sed$SppAlias == "Phy"),],
       aes(y=`Diet_up%`*100, x = C_s/C_WTO,
           group = interaction(loc, SppAlias, PFAA),
           color = PFAA))+
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed")+
  scale_color_manual(values = PFAA_cols)+
  facet_wrap(loc~SppAlias)+
  # scale_linetype_manual(values = c("dashed", "dotted", "solid"))+
  # ylim(-7,100)+
  # xlim(0,1)+
  geom_line(size= 0.5, alpha = 0.5)+
  theme_bw()+
  xlab("PFAS in sediment:water ratio (ng/g: ng/mL)")+
  ylab("PFAS dietary uptake (%)")

p_sed2
# ggsave(filename = "../Figures/fig5_sedRatio_oct31.png", p_sed2, width = 6, height = 5)


