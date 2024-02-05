
library(dplyr)
library(tidyr)
library(ggplot2)

# a function to run food web model
source("./Code/PFAS_model_Sun_etal_classes_KK.R") # Set up classes
source("./Code/PFAS_model_Sun_etal_steadyStateMod_KK.R")
source("./Code/PFAS_model_Sun_etal_BioaccumulationMod.R")
source("./Code/runEcosystemModels.R")
options(error = traceback)

# All runs go through model estimates without renal estimates 
# Set up function to run model with and without renal correction ----
# ********************************************************
# ********************************************************
GE_inputFiles_list = list(
    'numSpecies'= 12,
    'organismData'= read.csv('./Data/Sun_etal/organismData_Munoz.csv', row.names = "X"),
    'chemicalData'= read.csv('./Data/Sun_etal/chemicalData_Munoz.csv', row.names = "chemicalParameter"),
    'oceanData'= read.csv('./Data/Sun_etal/oceanData_Munoz.csv', row.names = "oceanParameter"),
    'chemicalParams'= read.csv('./Data/Sun_etal/chemicalParameters.csv', row.names = "chemicalParameter"),
    'foodWebData'= read.csv('./Data/Sun_etal/foodWebTable_Munoz.csv', row.names = "X")
)

GE_dietData = read.csv('./Data/Sun_etal/dietData_median_Munoz.csv', row.names = "X")  # changed to csv
GE_kRTable = read.csv('./Data/export_krTable.csv', row.names = "X", header = T)
colnames(GE_kRTable)<-c("kr/kb", "PFAA", "chemID")
GE_min_dietData = read.csv('./Data/Sun_etal/dietData_min_Munoz.csv', row.names = "X")
GE_max_dietData = read.csv('./Data/Sun_etal/dietData_max_Munoz.csv', row.names = "X")

# VERSION 1: food web - forced diets----
GE_settings_obsDiet = new_Settings(
  chooseStudyType = 'default',
  chooseDiet = 'forced Munoz',
  chooseEd = 'Goeritz',
  chooseKoc = 'Koc_Munoz',
  chooseRenal ='off'
)

AllData<-runEcosystemModel(settings = GE_settings_obsDiet,
                  inputFiles_list = GE_inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k2','ke','kg','Total Elimination'),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = GE_dietData,
                  kRTable = GE_kRTable,
                  run_w_RenalElim = TRUE,
                  food_web_calc = TRUE,
                  min_diet = TRUE,
                  max_diet = TRUE,
                  min_dietData = GE_min_dietData, 
                  max_dietData = GE_max_dietData)

# Remove data for cases in which > 50% are NDs
# instead of deleting these rows, replacing with np.nan allows them to be skipped over in plots
AllData[c(c(AllData$SppAlias=='Fln' & c(AllData$PFAA == 'PFOA' | AllData$PFAA == 'PFHxS')) |
                   c(AllData$SppAlias=='Gob' & c(AllData$PFAA == 'PFOA' | AllData$PFAA =='PFHxS')) |
                   c(AllData$SppAlias=='CSb' & AllData$PFAA=='PFHxS') |
                   c(AllData$SppAlias=='Acy' & AllData$PFAA=='PFOA')|
                   c(AllData$SppAlias=='Spr')), c(3:ncol(AllData))] <- NaN


# VERSION 2: full food web ----
GE_settings_modelDiet = new_Settings(
  chooseStudyType = 'default',
  chooseDiet = 'default',
  chooseEd = 'Goeritz',
  chooseKoc = 'Koc_Munoz',
  chooseRenal ='off'
)

AllData_fw<-runEcosystemModel(settings = GE_settings_modelDiet,
                  inputFiles_list = GE_inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k2','ke','kg','Total Elimination'),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = GE_dietData,
                  kRTable = GE_kRTable,
                  run_w_RenalElim = TRUE,
                  food_web_calc = TRUE,
                  min_diet = TRUE,
                  max_diet = TRUE,
                  min_dietData = GE_min_dietData, 
                  max_dietData = GE_max_dietData)

# Remove data for cases in which > 50% are NDs
AllData_fw[c(c(AllData$SppAlias=='Fln' & c(AllData$PFAA == 'PFOA' | AllData$PFAA == 'PFHxS')) |
                   c(AllData$SppAlias=='Gob' & c(AllData$PFAA == 'PFOA' | AllData$PFAA =='PFHxS')) |
                   c(AllData$SppAlias=='CSb' & AllData$PFAA=='PFHxS') |
                   c(AllData$SppAlias=='Acy' & AllData$PFAA=='PFOA') |
                   c(AllData$SppAlias=='Spr')), c(3:ncol(AllData))] <- NaN


# ********************************************************
# Figure: Food web --------
# A. Plot MUNOZ CONCENTRATIONS 1:1
# Select species of interest
PlotData <- AllData
PlotData_fw <- AllData_fw

# Subset Data
FishData <- PlotData[c(PlotData$SppAlias %in% c('Sol','Fln','Gob','CSb','Acy','Spr')), ]
InvData <- PlotData[c(PlotData$SppAlias %in% c('Cop','Mys','Rag','WSh','Gam')), ]
FishData_fw <- PlotData_fw[c(PlotData$SppAlias %in% c('Sol','Fln','Gob','CSb','Acy','Spr')), ]
InvData_fw <- PlotData_fw[c(PlotData$SppAlias %in% c('Cop','Mys','Rag','WSh','Gam')), ]

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

p1_inv<-ggplot(data = InvData, aes(x = log_ngkg, y = Obs_logngkg, color = PFAA))+
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
    geom_errorbarh(aes(xmin = log_ngkg - logngkg_Lo, xmax = log_ngkg + logngkg_Up), size = 0.5) +
    geom_point(pch = 19, size = 2)+
    theme_classic()+
    scale_color_manual(values = PFAA_cols)+
    ggtitle(label = 'Inverts, Empirical diet data')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    coord_cartesian(ylim = c(0,5.5), xlim = c(0,5.5))+
    theme(legend.position = "none")+
    theme(legend.position = "none")

p1_inv_re<-ggplot(data = InvData_fw, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
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
    ggtitle(label = 'Inverts, Simulated diet data')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    coord_cartesian(ylim = c(0,5.5), xlim = c(0,5.5))+
    theme(legend.position = "none")

p1_fish<-ggplot(data = FishData, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
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
    ggtitle(label = 'Fish, Empirical diet data!!??')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    coord_cartesian(ylim = c(0,5.5), xlim = c(0,5.5))+
    theme(legend.position = "none")
# p1_fish

p1_fish_re<-ggplot(data = FishData_fw, aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
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
    coord_cartesian(ylim = c(0,5.5), xlim = c(0,5.5))+
    # theme(legend.position = "right")
    theme(legend.position = "none")
# p1_fish_re
legend_plot <- ggplot(data = FishData_fw,
         aes(x = log_ngkg_re, y = Obs_logngkg, color = PFAA))+
  geom_point()+
  scale_color_manual(values = PFAA_cols)+
  theme(legend.title = element_blank())+
  theme_classic()

legend<-cowplot::get_legend(legend_plot)
fig3<-cowplot::plot_grid(p1_inv, p1_inv_re, legend,
                   p1_fish, p1_fish_re, NULL,
                   ncol = 3,
                   nrow = 2,
                   rel_widths = c(1,1, 0.3),
        labels = c("A", "B", "", "C", "D", ""))
# ggsave(filename = "./Figures/Figure3_inv_fish.png", plot = fig3,  width = 7, height = 6.5)




