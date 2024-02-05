
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

# FALL SEASON -----
# All runs go through model estimates without renal estimates 
## 1. Set up function to run model with and without renal correction ----
# ********************************************************
# ********************************************************
WG_inputFiles_list = list(
    'numSpecies'= 4,
    'organismData'= read.csv('./Data/WG/organismData_BMFBCF_WG_fall.csv', row.names = "X"),
    'chemicalData'= read.csv('./Data/WG/chemicalData_full_WG_fall.csv', row.names = "chemicalParameter"),
    'oceanData'= read.csv('./Data/WG/oceanData_BCFBMF_WG_fall.csv', row.names = "oceanParameter"),
    'chemicalParams'= read.csv('./Data/WG/chemicalParameters_BCFBMF_WG.csv', row.names = "chemicalParameter"),
    'foodWebData'= read.csv('./Data/WG/foodWebTable_BCFBMF_WG_fish_fall.csv', row.names = "X")
)

WG_kRTable = read.csv('./Data/export_krTable.csv', row.names = "X", header = T)
colnames(WG_kRTable)<-c("kr/kb", "PFAA", "chemID")

WG_dietData = read.csv('./Data/WG/dietData_median_WG_fall.csv', row.names = "Species")  # changed to csv
WG_min_dietData = read.csv('./Data/WG/dietData_min_WG_fall.csv', row.names = "Species")
WG_max_dietData = read.csv('./Data/WG/dietData_max_WG_fall.csv', row.names = "Species")

## 2. food web - forced diets----
WG_settings_obsDiet = new_Settings(
  chooseStudyType = 'default',
  chooseDiet = 'forced Munoz',
  chooseEd = 'Goeritz',
  chooseKoc = 'Koc_Munoz',
  chooseRenal ='off'
)

AllData<-runEcosystemModel(settings = WG_settings_obsDiet,
                  inputFiles_list = WG_inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k2','ke','kg','Total Elimination'),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = WG_dietData,
                  kRTable = WG_kRTable,
                  run_w_RenalElim = TRUE,
                  food_web_calc = TRUE,
                  min_diet = TRUE,
                  max_diet = TRUE,
                  min_dietData = WG_min_dietData, 
                  max_dietData = WG_max_dietData)



## 3. food web - predicted ----
WG_settings_modelDiet = new_Settings(
  chooseStudyType = 'default',
  chooseDiet = 'default',
  chooseEd = 'Goeritz',
  chooseKoc = 'Koc_Munoz',
  chooseRenal ='off'
)

AllData_fw<-runEcosystemModel(settings = WG_settings_modelDiet,
                  inputFiles_list = WG_inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k2','ke','kg','Total Elimination'),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = WG_dietData,
                  kRTable = WG_kRTable,
                  run_w_RenalElim = TRUE,
                  food_web_calc = TRUE,
                  min_diet = TRUE,
                  max_diet = TRUE,
                  min_dietData = WG_min_dietData, 
                  max_dietData = WG_max_dietData)




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
fig3<-cowplot::plot_grid(p1, p1_re, legend,
                   ncol = 3,
                   nrow = 1,
                   rel_widths = c(1,1, 0.3),
        labels = c("A", "B", ""))
ggsave(filename = "./Figures/Figure3_WG_fish_FALL.png", plot = fig3,  width = 10, height = 4)
fig3



# ********************************************************




# Other Figures --------

ggplot(AllData_fw, aes(PFAA, log(Cb_change), color = PFAA, fill= PFAA))+
  geom_boxplot(alpha=0.1)+
  geom_point(size=3, pch=21, color = "black")+
  scale_color_manual(values = PFAA_cols)+
  scale_fill_manual(values = PFAA_cols)+
  theme_classic()+
  ggtitle("Change in log(tissue concentration) w/ and wo/ Renal: \nfall")+
  theme(legend.position = "none")

data<-cbind(AllData_fw[, c("PFAA", "log_ngkg_re", "logngkg_re_Lo","logngkg_re_Up", "SppAlias")],
     AllData[, c("PFAA", "log_ngkg_re", "logngkg_re_Lo","logngkg_re_Up", "SppAlias")])

colnames(data)<-c("PFAA", "log_ngkg_re", "logngkg_re_Lo","logngkg_re_Up", "SppAlias",
                          "PFAA_emp", "log_ngkg_re_emp", "logngkg_re_Lo_emp","logngkg_re_Up_emp", "SppAlias_emp")


ggplot(data, aes(log_ngkg_re, log_ngkg_re_emp, color = PFAA, label = SppAlias))+
  geom_point()+
  scale_color_manual(values = PFAA_cols)+
  geom_text(check_overlap = T)+
  coord_cartesian(ylim = c(2,7), xlim = c(2,7))

ggplot(data, aes(logngkg_re_Up, logngkg_re_Up_emp, color = PFAA, label = SppAlias))+
  geom_point()+
  scale_color_manual(values = PFAA_cols)+
  geom_text(check_overlap = T)+
  coord_cartesian(ylim = c(0,0.7), xlim = c(0,0.7))

ggplot(data, aes(logngkg_re_Lo, logngkg_re_Lo_emp, color = PFAA, label = SppAlias))+
  geom_point()+
  scale_color_manual(values = PFAA_cols)+
  geom_text(check_overlap = T)+
  coord_cartesian(ylim = c(0,5), xlim = c(0,5))



  




