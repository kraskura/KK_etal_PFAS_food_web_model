
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)

# a function to run food web model
source("./Code/PFAS_model_Sun_etal_classes_KK.R") # Set up classes
source("./Code/PFAS_model_Sun_etal_steadyStateMod_KK.R")
source("./Code/PFAS_model_Sun_etal_BioaccumulationMod.R")
source("./Code/runEcosystemModels.R")
options(error = traceback)

# SPRING SEASON -----
# All runs go through model estimates without renal estimates 
## 1. Set up function to run model with and without renal correction ----
# ********************************************************
# ********************************************************
JBA_inputFiles_list = list(
    'numSpecies'= 9,
    'organismData'= read.csv('./Data/JBA/organismData_BMFBCF_JBA_spring.csv', row.names = "X"),
    'chemicalData'= read.csv('./Data/JBA/chemicalData_full_JBA_spring.csv', row.names = "chemicalParameter"),
    'oceanData'= read.csv('./Data/JBA/oceanData_BCFBMF_JBA_spring.csv', row.names = "oceanParameter"),
    'chemicalParams'= read.csv('./Data/JBA/chemicalParameters_BCFBMF_JBA.csv', row.names = "chemicalParameter"),
    'foodWebData'= read.csv('./Data/JBA/foodWebTable_BCFBMF_JBA_fish_springV2.csv', row.names = "X")
)

JBA_kRTable = read.csv('./Data/export_krTable.csv', row.names = "X", header = T)
colnames(JBA_kRTable)<-c("kr/kb", "PFAA", "chemID")

JBA_dietData = read.csv('./Data/JBA/dietData_median_JBA_spring.csv', row.names = "Species")  # changed to csv
JBA_min_dietData = read.csv('./Data/JBA/dietData_min_JBA_spring.csv', row.names = "Species")
JBA_max_dietData = read.csv('./Data/JBA/dietData_max_JBA_spring.csv', row.names = "Species")

## 2. food web - forced diets----
JBA_settings_obsDiet = new_Settings(
  chooseStudyType = 'default',
  chooseDiet = 'forced Munoz',
  chooseEd = 'Goeritz',
  chooseKoc = 'Koc_Munoz',
  chooseRenal ='off'
)

AllData<-runEcosystemModel(settings = JBA_settings_obsDiet,
                  inputFiles_list = JBA_inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k2','ke','kg','Total Elimination'),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = JBA_dietData,
                  kRTable = JBA_kRTable,
                  run_w_RenalElim = TRUE,
                  food_web_calc = TRUE,
                  min_diet = TRUE,
                  max_diet = TRUE,
                  min_dietData = JBA_min_dietData, 
                  max_dietData = JBA_max_dietData)



## 3. food web - predicted ----
JBA_settings_modelDiet = new_Settings(
  chooseStudyType = 'default',
  chooseDiet = 'default',
  chooseEd = 'Goeritz',
  chooseKoc = 'Koc_Munoz',
  chooseRenal ='off'
)

AllData_fw<-runEcosystemModel(settings = JBA_settings_modelDiet,
                  inputFiles_list = JBA_inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k2','ke','kg','Total Elimination'),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = JBA_dietData,
                  kRTable = JBA_kRTable,
                  run_w_RenalElim = TRUE,
                  food_web_calc = TRUE,
                  min_diet = TRUE,
                  max_diet = TRUE,
                  min_dietData = JBA_min_dietData, 
                  max_dietData = JBA_max_dietData)




# ********************************************************
## Figure: Food web --------
# A. Plot MUNOZ CONCENTRATIONS 1:1
# Select species of interest
PlotData <- AllData
PlotData_fw <- AllData_fw

# Set up PFAAs and labels for plot looping
PFAS_List <- c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA')
# PFAS_labels <- c('C6 PFSA','C8 PFSA','C8 PFCA', 'C9 PFCA','C10 PFCA','C11 PFCA')
PFAA_cols <- c('PFHxS' = 'blueviolet', # 'mediumslateblue' not plotted?
             'PFOS' = 'cornflowerblue', # 'c'
             'PFOA' = 'orange',
             'PFNA' = 'violet',
             'PFDA' = 'deeppink',
             'PFUA' = 'limegreen')

p_lines <- data.frame('x' = c(0, 7),
                      'y' = c(0, 7))
p_lines$y1 <- p_lines$x-1
p_lines$y2 <- p_lines$x+1
p_lines$y3 <- p_lines$x-log10(2)
p_lines$y4 <- p_lines$x+log10(2)

# DataList$PFAA<-factor(DataList$PFAA)

p1<-ggplot(data = AllData, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
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
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_Lo, xmax = log_ngkg_re + logngkg_Up), size = 0.5) +
    geom_point(pch = 19, size = 2)+
    theme_classic()+
    facet_grid(.~SppAlias)+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Fish, empirical data')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(2,6.5), xlim = c(2,6.5))+
    theme(legend.position = "none")
# p1_fish

p1_re<-ggplot(data = AllData_fw, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
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
    geom_point(pch = 19, size = 2)+
    theme_classic()+
    facet_grid(.~SppAlias)+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Fish, Simulated diet data')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(0,6.5), xlim = c(0,6.5))+
    # theme(legend.position = "right")
    theme(legend.position = "none")
# p1_fish_re
legend_plot <- ggplot(data = AllData_fw,
         aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
  geom_point()+
  scale_color_manual(values = PFAA_cols)+
  theme(legend.title = element_blank())+
  theme_classic()

legend<-cowplot::get_legend(legend_plot)
fig3<-cowplot::plot_grid(p1,  legend,p1_re,
                   ncol = 2,
                   nrow = 2,
                   rel_widths = c(1.9, 0.3),
        labels = c("A", "B", ""))
ggsave(filename = "./Figures/Figure3_JBA_fish_SPRING_foodwebV2.png", plot = fig3,  width = 14, height = 5)


# SUMMER SEASON -----

# All runs go through model estimates without renal estimates 
## 1. Set up function to run model with and without renal correction ----
# ********************************************************
# ********************************************************
JBA_inputFiles_list_summer = list(
    'numSpecies'= 9,
    'organismData'= read.csv('./Data/JBA/organismData_BMFBCF_JBA_summer.csv', row.names = "X"),
    'chemicalData'= read.csv('./Data/JBA/chemicalData_full_JBA_summer.csv', row.names = "chemicalParameter"),
    'oceanData'= read.csv('./Data/JBA/oceanData_BCFBMF_JBA_summer.csv', row.names = "oceanParameter"),
    'chemicalParams'= read.csv('./Data/JBA/chemicalParameters_BCFBMF_JBA.csv', row.names = "chemicalParameter"),
    'foodWebData'= read.csv('./Data/JBA/foodWebTable_BCFBMF_JBA_fish_summerV2.csv', row.names = "X")
)

JBA_kRTable = read.csv('./Data/export_krTable.csv', row.names = "X", header = T)
colnames(JBA_kRTable)<-c("kr/kb", "PFAA", "chemID")

JBA_dietData_summer = read.csv('./Data/JBA/dietData_median_JBA_summer.csv', row.names = "Species")  # changed to csv
JBA_min_dietData_summer = read.csv('./Data/JBA/dietData_min_JBA_summer.csv', row.names = "Species")
JBA_max_dietData_summer = read.csv('./Data/JBA/dietData_max_JBA_summer.csv', row.names = "Species")

## 2. food web - forced diets----
JBA_settings_obsDiet_summer = new_Settings(
  chooseStudyType = 'default',
  chooseDiet = 'forced Munoz',
  chooseEd = 'Goeritz',
  chooseKoc = 'Koc_Munoz',
  chooseRenal ='off'
)

AllData_summer<-runEcosystemModel(settings = JBA_settings_obsDiet_summer,
                  inputFiles_list = JBA_inputFiles_list_summer,
                  parameterList = c('C_B','C_WTO','C_s','k2','ke','kg','Total Elimination'),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = JBA_dietData_summer,
                  kRTable = JBA_kRTable,
                  run_w_RenalElim = TRUE,
                  food_web_calc = TRUE,
                  min_diet = TRUE,
                  max_diet = TRUE,
                  min_dietData = JBA_min_dietData_summer, 
                  max_dietData = JBA_max_dietData_summer)



## 3. food web - predicted ----
JBA_settings_modelDiet_summer = new_Settings(
  chooseStudyType = 'default',
  chooseDiet = 'default',
  chooseEd = 'Goeritz',
  chooseKoc = 'Koc_Munoz',
  chooseRenal ='off'
)

AllData_fw_summer<-runEcosystemModel(settings = JBA_settings_modelDiet_summer,
                  inputFiles_list = JBA_inputFiles_list_summer,
                  parameterList = c('C_B','C_WTO','C_s','k2','ke','kg','Total Elimination'),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = JBA_dietData_summer,
                  kRTable = JBA_kRTable,
                  run_w_RenalElim = TRUE,
                  food_web_calc = TRUE,
                  min_diet = TRUE,
                  max_diet = TRUE,
                  min_dietData = JBA_min_dietData_summer, 
                  max_dietData = JBA_max_dietData_summer)




# ********************************************************
## Figure: Food web --------

p1_summer<-ggplot(data = AllData_summer, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
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
    geom_errorbarh(aes(xmin = log_ngkg_re - logngkg_Lo, xmax = log_ngkg_re + logngkg_Up), size = 0.5) +
    geom_point(pch = 19, size = 2)+
    theme_classic()+
    facet_grid(.~SppAlias)+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Fish, empirical data')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(2,6.5), xlim = c(2,6.5))+
    theme(legend.position = "none")
# p1_fish

p1_re_summer<-ggplot(data = AllData_fw_summer, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
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
    geom_point(pch = 19, size = 2)+
    theme_classic()+
    facet_grid(.~SppAlias)+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Fish, Simulated diet data')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    # coord_cartesian(ylim = c(0,6.5), xlim = c(0,6.5))+
    # theme(legend.position = "right")
    theme(legend.position = "none")
# p1_fish_re

# p1_summerNOREN<-ggplot(data = AllData_fw_summer, aes(x = log_ngkg, y = Obs_logngkg, color = PFAA))+
#     geom_line(data = p_lines, mapping = aes(x = x, y = y),
#               linetype = 1, linewidth = 0.5, color = "black")+
#     geom_line(data = p_lines, mapping = aes(x = x, y = y1),
#               linetype = 3, linewidth = 0.2, color = "black")+
#     geom_line(data = p_lines, mapping = aes(x = x, y = y2),
#               linetype = 3, linewidth = 0.2, color = "black")+
#     geom_line(data = p_lines, mapping = aes(x = x, y = y3),
#               linetype = 3, linewidth = 0.2, color = "black")+
#     geom_line(data = p_lines, mapping = aes(x = x, y = y4),
#               linetype = 3, linewidth = 0.2, color = "black")+
#     geom_errorbar(aes(ymin = Obs_logngkg - Obs_logngkg_Lo, ymax = Obs_logngkg + Obs_logngkg_Up), size = 0.5) +
#     geom_errorbarh(aes(xmin = log_ngkg - logngkg_Lo, xmax = log_ngkg + logngkg_Up), size = 0.5) +
#     geom_point(pch = 19, size = 2)+
#     theme_classic()+
#     scale_color_manual(values = PFAA_cols)+
#     ggtitle(label = 'Fish, Simulated diet data')+
#     xlab('Modeled log PFAA (ng/kg)')+
#     ylab('Observed log PFAA (ng/kg)')+
#     # coord_cartesian(ylim = c(0,6.5), xlim = c(0,6.5))+
#     # theme(legend.position = "right")
#     theme(legend.position = "none")

legend<-cowplot::get_legend(legend_plot)

fig3_summer<-cowplot::plot_grid(p1_summer,  legend, p1_re_summer,
                   ncol = 2,
                   nrow = 2,
                   rel_widths = c(1.9, 0.3),
        labels = c("A", "B", ""))
ggsave(filename = "./Figures/Figure3_JBA_fish_SUMMER_foodwebV2.png", plot = fig3_summer,  width = 14, height = 5)
fig3_summer



# Other Figures --------
ggplot(AllData_fw_summer, aes(PFAA, log(Cb_change), color = PFAA, fill= PFAA))+
  geom_boxplot(alpha=0.1)+
  geom_point(size=3, pch=21, color = "black")+
  scale_color_manual(values = PFAA_cols)+
  scale_fill_manual(values = PFAA_cols)+
  theme_classic()+
  ggtitle("Change in log(tissue concentration) w/ and wo/ Renal: \nSUMMER")+
  theme(legend.position = "none")

ggplot(AllData_fw, aes(PFAA, log(Cb_change), color = PFAA, fill= PFAA))+
  geom_boxplot(alpha=0.1)+
  geom_point(size=3, pch=21, color = "black")+
  scale_color_manual(values = PFAA_cols)+
  scale_fill_manual(values = PFAA_cols)+
  theme_classic()+
  ggtitle("Change in log(tissue concentration) w/ and wo/ Renal: \nSPRING")+
  theme(legend.position = "none")

data_summer<-cbind(AllData_fw_summer[, c("PFAA", "log_ngkg_re", "logngkg_re_Lo","logngkg_re_Up","Obs_logngkg",  "SppAlias")],
     AllData_summer[, c("PFAA", "log_ngkg_re", "logngkg_re_Lo","logngkg_re_Up","SppAlias")])

colnames(data_summer)<-c("PFAA", "log_ngkg_re", "logngkg_re_Lo","logngkg_re_Up", "Obs_logngkg","SppAlias",
                          "PFAA_emp", "log_ngkg_re_emp", "logngkg_re_Lo_emp","logngkg_re_Up_emp", "SppAlias_emp")

data_spring<-cbind(AllData_fw[, c("PFAA", "log_ngkg_re", "logngkg_re_Lo","logngkg_re_Up","Obs_logngkg", "SppAlias")],
     AllData[, c("PFAA", "log_ngkg_re", "logngkg_re_Lo","logngkg_re_Up", "SppAlias")])

colnames(data_spring)<-c("PFAA", "log_ngkg_re", "logngkg_re_Lo","logngkg_re_Up","Obs_logngkg", "SppAlias",
                          "PFAA_emp", "log_ngkg_re_emp", "logngkg_re_Lo_emp","logngkg_re_Up_emp", "SppAlias_emp")


ggplot(data_summer, aes(log_ngkg_re, Obs_logngkg, color = PFAA, label = SppAlias))+
  geom_point()+
  scale_color_manual(values = PFAA_cols)+
  geom_text(check_overlap = T)+
  coord_cartesian(ylim = c(2,7), xlim = c(2,7))

ggplot(data_spring, aes(log_ngkg_re, Obs_logngkg, color = PFAA, label = SppAlias))+
  geom_point()+
  scale_color_manual(values = PFAA_cols)+
  geom_text(check_overlap = T)+
  coord_cartesian(ylim = c(2,7), xlim = c(2,7))

ggplot(data_summer, aes(log_ngkg_re, log_ngkg_re_emp, color = PFAA, label = SppAlias))+
  geom_point()+
  scale_color_manual(values = PFAA_cols)+
  geom_text(check_overlap = T)+
  coord_cartesian(ylim = c(0,5), xlim = c(0,5))



  




