# date: march 28 2024
# author:  Krista Kraskura
# description: organizing MDE data

# setup, libraries ---- 
library(tidyverse)
library(here)
library(readxl)
library(lubridate)

# PFAA acronyms 



# data wrangling -----
# read in their data 
d_fish_21<-read_excel(here("Data/MDE/2021AdvisoryData_PFAS_AWQMSFinal.xlsx"), sheet = 1)
d_fish_22<-read_excel(here("Data/MDE/2022AdvisoryData_PFAS_AWQMSFinal.xlsx"), sheet = 1)
d_fish_23<-read_excel(here("Data/MDE/2023AdvisoryData_GEL_PFAS_AWQMSFinal.xlsx"), sheet = 1)

d_fishPisc_22<-read_excel(here("Data/MDE/2022AdvisoryData_PFAS_AWQMSFinal_Piscataway.xlsx"), sheet = 1)
d_fish<-read_excel(here("Data/MDE/PFAS_Tissue_2020_21_22.xlsx"), sheet = 1, skip = 1)

d_water<-read_excel(here("Data/MDE/PFAS_SurfaceWater_2021_22.xlsx"), sheet = 1, skip = 1)

# are data in d_fish 21, 22, and 23 overlap with data in Piscataway data and PFAS_Tissue_2020_21_22 data?
# identify duplicates 
n21<-length(unique(d_fish_21$Activity_ID)) # 19
n22<-length(unique(d_fish_22$Activity_ID)) # 68
n23<-length(unique(d_fish_23$Activity_ID)) # 65
n20s<-length(unique(d_fish$`Activity ID`)) # 152
np<-length(unique(d_fishPisc_22$Activity_ID)) # 8

n21+n22+n23+n20s+np # 312 
n21+n22+n23 == n20s # this is a combo of all years?

all(c((levels(factor(d_fish_21$Activity_ID))),
  (levels(factor(d_fish_22$Activity_ID))), 
  (levels(factor(d_fish_23$Activity_ID))))  ==
  levels(factor(d_fish$`Activity ID`))) # all Activity IDs are the same 

# do Piscataway fish are also in the 20-22 data? --> no
grepl(pattern = paste(levels(factor(d_fishPisc_22$Activity_ID)), collapse = "|"),
      x = levels(factor(d_fish$`Activity ID`)), ignore.case = TRUE)

# merge piscataway and other fish data
d_fish<-rbind(d_fish_21, d_fish_22, d_fish_23, d_fishPisc_22)

# format all sampling dates as dates
d_fish$ActivitySdate<-as.Date(d_fish$ActivitySdate)
d_water$`Activity Start Date`<-as.Date(d_water$`Activity Start Date`)

# take out all columns that contain NA only
d_fish<-d_fish[,which(unlist(lapply(d_fish, function(x)!all(is.na(x))))),with=F]
d_water<-d_water[,which(unlist(lapply(d_water, function(x)!all(is.na(x))))),with=F]

# get matching column names to combine all data in one sheet
names(d_fish) <- gsub(" ", "_", names(d_fish))
names(d_water) <- gsub(" ", "_", names(d_water))

# need project ID
d_fish$Organization_ID<-"MDE"
names(d_fish)[names(d_fish) == "Monitoring_L_ID"] <- "Monitoring_Location_ID"
names(d_fish)[names(d_fish) == "ActivitySdate"] <- "Activity_Start_Date"
names(d_fish)[names(d_fish) == "Project_name"] <- "Project_ID"
names(d_fish)[names(d_fish) == "RD_Q_limit_measure"] <- "Detection_Limit_Value"
names(d_fish)[names(d_fish) == "RD_Q_limit_unit"] <- "Detection_Limit_Unit"
names(d_fish)[names(d_fish) == "Lab_name"] <- "Laboratory_Name"
names(d_fish)[names(d_fish) == "Activity_Media_Name"] <- "Activity_Media"
names(d_fish)[names(d_fish) == "Char_name"] <- "Characteristic_Name"
names(d_water)[names(d_water) == "Detection_Limit_Value1"] <- "Detection_Limit_Value"
names(d_water)[names(d_water) == "Detection_Limit_Unit1"] <- "Detection_Limit_Unit"
names(d_water)[names(d_water) == "Project_ID1"] <- "Project_ID"
d_water$Subject_Taxonomic_Name <- "water" # no species in this data

data<-rbind(d_fish[ ,c("Organization_ID","Project_ID", "Activity_ID", "Monitoring_Location_ID",
                "Activity_Type", "Activity_Start_Date", "Laboratory_Name", "Result_Value_Type", "Result_Value",
                "Subject_Taxonomic_Name","Characteristic_Name", "Activity_Media",
                "Sample_Collection_Method_ID", "Detection_Limit_Unit", "Detection_Limit_Value")],
      d_water[ ,c("Organization_ID","Project_ID", "Activity_ID", "Monitoring_Location_ID",
                "Activity_Type","Activity_Start_Date", "Laboratory_Name", "Result_Value_Type", "Result_Value",
                "Subject_Taxonomic_Name","Characteristic_Name", "Activity_Media",
                "Sample_Collection_Method_ID", "Detection_Limit_Unit", "Detection_Limit_Value")])

# make all lower case
data$Characteristic_Name<-tolower(data$Characteristic_Name)
levels(factor(data$Characteristic_Name))

# subset data 
data_m<-data %>% 
  filter(Characteristic_Name == "average weight") # body size 
data_len<-data %>% 
  filter(Characteristic_Name == "average length") # body length
data_n<-data %>% 
  filter(Characteristic_Name == "number of individuals") # number of individuals
data_tot<-data %>% 
  filter(Characteristic_Name == "total pfas") # number of individuals

data_l<-data %>%  # pfas data
  filter(!Characteristic_Name == "average weight", 
         !Characteristic_Name == "average length", 
         !Characteristic_Name == "number of individuals", 
         !Characteristic_Name == "total pfas")

# n of unique sample sizes
ntot<-length(unique(data_l$Activity_ID)) # 242 samples 
ntot_m<-length(unique(data_m$Activity_ID)) # 159
ntot_l<-length(unique(data_l$Activity_ID)) # 242
ntot_n<-length(unique(data_n$Activity_ID)) # 160

names(data_m)[names(data_m) == 'Result_Value'] <- 'Mass_g'
names(data_len)[names(data_len) == 'Result_Value'] <- 'Length_cm'
names(data_n)[names(data_n) == 'Result_Value'] <- 'n'
names(data_tot)[names(data_tot) == 'Result_Value'] <- 'Sum_PFAS'

# merge data sets and format
data_l<-(merge(data_l, data_m[, c("Activity_ID", "Mass_g")], by = "Activity_ID", all.x = TRUE))
data_l<-(merge(data_l, data_len[, c("Activity_ID", "Length_cm")], by = "Activity_ID", all.x = TRUE))
data_l<-(merge(data_l, data_n[, c("Activity_ID", "n")], by = "Activity_ID", all.x = TRUE))
data_l<-(merge(data_l, data_tot[, c("Activity_ID", "Sum_PFAS")], by = "Activity_ID", all.x = TRUE))

data_l<-data_l %>% 
  mutate_at(c("Result_Value", "Detection_Limit_Value", "Length_cm", "Mass_g", "n", "Sum_PFAS"), as.numeric) %>% 
  mutate_at(c("Activity_ID", "Organization_ID", "Project_ID", "Monitoring_Location_ID", 
              "Activity_Type", "Laboratory_Name", "Result_Value_Type", 
              "Subject_Taxonomic_Name", "Characteristic_Name", "Activity_Media", 
              "Sample_Collection_Method_ID", "Detection_Limit_Unit"), as.factor) 

# percent PFAS
data_l$pct_PFAS<-data_l$Result_Value/data_l$Sum_PFAS
data_l$Year<-lubridate::year(data_l$Activity_Start_Date)
data_l$Month<-lubridate::month(data_l$Activity_Start_Date)
# pivot wide format data
data_w<-data_l %>% 
  pivot_wider(names_from = Characteristic_Name,
              values_from = Result_Value)


# Figures ------ 

data_l_sum<-data_l %>% 
  dplyr::group_by(Year, Monitoring_Location_ID,
                  Characteristic_Name) %>% 
  summarise(Result_Value_mean = mean(Result_Value),
            pct_PFAS_mean = mean(pct_PFAS),
            n_mean = n()) 
  

year <- 2022


  data_l.w <- data_l %>%
    filter(Year == year & 
           Subject_Taxonomic_Name == "water")
  levels(data_l.w$Subject_Taxonomic_Name)
  levels(data_l.w$Monitoring_Location_ID)

  data_l.f <- data_l %>%
    filter(Year == year & 
         !Subject_Taxonomic_Name == "water")
  levels(data_l$Subject_Taxonomic_Name)
  levels(data_l.w$Monitoring_Location_ID)
    
  p1<-ggplot(data_l.f,
             aes(Monitoring_Location_ID, Result_Value,
                  group = Monitoring_Location_ID,
                  color = Characteristic_Name,
                  fill = Characteristic_Name))+
    geom_bar(stat = "identity",  show.legend = F)+
    theme_bw()+
    scale_fill_viridis_d()+
    scale_color_viridis_d()+
    theme(axis.text.x = element_text(angle = 90, face = "italic"))+
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          panel.spacing = unit(0, "lines"))+
    xlab("Species in each sampling area")+
    ylab("PFAS (ng/g)")
  
  
  p1.w<-ggplot(data_l.w,
               aes(Monitoring_Location_ID, Result_Value,
                    group = Monitoring_Location_ID,
                    color = Characteristic_Name,
                    fill = Characteristic_Name))+
      geom_bar(stat = "identity",  show.legend = F)+
      # facet_grid(.~Monitoring_Location_ID, switch="both")+
      theme_bw()+
      scale_fill_viridis_d()+
      scale_color_viridis_d()+
      theme(axis.text.x = element_text(angle = 90, face = "italic"))+
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            panel.spacing = unit(0, "lines"))+
      xlab("Species in each sampling area")+
      ylab("PFAS (ng/g)")+
      # ylim(-250, 250)+
      geom_bar(data_l.f, 
               mapping = aes(Monitoring_Location_ID, -Result_Value,
                    group = Monitoring_Location_ID,
                    color = Characteristic_Name,
                    fill = Characteristic_Name),
        stat = "identity",  show.legend = F)
  
  # p2<-ggplot(data_l_sum[data_l_sum$Year == year, ],
  #            aes(Monitoring_Location_ID, pct_PFAS_mean,
  #                 group = Monitoring_Location_ID,
  #                 color = Characteristic_Name,
  #                 fill = Characteristic_Name))+
  #   geom_bar(stat = "identity", show.legend = F) +
  #   # facet_grid(.~Monitoring_Location_ID, switch="both")+
  #   theme_bw()+
  #   scale_fill_viridis_d()+
  #   scale_color_viridis_d()+
  #   theme(axis.text.x = element_text(angle = 90, face = "italic"))+
  #   theme(strip.background = element_blank(),
  #         strip.placement = "outside",
  #         panel.spacing = unit(0, "lines"))+
  #   xlab("Species in each sampling area")+
  #   ylab("PFAS (ng/g)")
  
  cowplot::plot_grid(p1, p1.w, nrow = 2, align = "hv")
  
  
  p1.b<-ggplot(data_l[data_l$Year == year & !data_l$Characteristic_Name == "perfluorooctane sulfonic acid", ],
               aes(Subject_Taxonomic_Name, Result_Value,
                                            group = Monitoring_Location_ID,
                                            color = Characteristic_Name,
                                            fill = Characteristic_Name))+
    geom_bar(stat = "identity",  show.legend = F)+
    facet_grid(.~Monitoring_Location_ID, switch="both")+
    theme_bw()+
    scale_fill_viridis_d()+
    scale_color_viridis_d()+
    theme(axis.text.x = element_text(angle = 90, face = "italic"))+
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          panel.spacing = unit(0, "lines"))
  # 
  # p2.b<-ggplot(data_l_sum[data_l_sum$Year == year & !data_l_sum$Characteristic_Name == "perfluorooctane sulfonic acid", ],
  #              aes(Subject_Taxonomic_Name, pct_PFAS_mean,
  #                                           group = Monitoring_Location_ID,
  #                                           color = Characteristic_Name,
  #                                           fill = Characteristic_Name))+
  #   geom_bar(stat = "identity", show.legend = F) +
  #   facet_grid(.~Monitoring_Location_ID, switch="both")+
  #   theme_bw()+
  #   scale_fill_viridis_d()+
  #   scale_color_viridis_d()+
  #   theme(axis.text.x = element_text(angle = 90, face = "italic"))+
  #   theme(strip.background = element_blank(),
  #         strip.placement = "outside",
  #         panel.spacing = unit(0, "lines"))
  # 
  # 
  cowplot::plot_grid(p1.b, nrow = 2)
  

  
  
  
  
# All years all locations 

  ggplot(data_l, aes(x = Activity_Start_Date, y = Result_Value,
                                            color = Characteristic_Name,
                                            fill = Characteristic_Name))+
    geom_point(pch = 1, color = "black", stroke = 0.2, stat = "identity", show.legend = F) +
    geom_point(data= data_l[data_l$Subject_Taxonomic_Name == "water", ], 
               aes(x = Activity_Start_Date, y = Result_Value),
               pch = 1, color = "dodgerblue3", stroke = 0.5, size = 1, stat = "identity", show.legend = F) +
    geom_point(data= data_l[data_l$Characteristic_Name == "perfluorooctane sulfonic acid", ], 
               aes(x = Activity_Start_Date, y = Result_Value),
               pch = 19, color = "red3",size = 0.5, stat = "identity", show.legend = F) +
    theme_bw()+
    scale_fill_viridis_d()+
    scale_color_viridis_d()+
    theme(axis.text.x = element_text(angle = 90, face = "italic"))+
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          panel.spacing = unit(0, "lines"))
  
    ggplot(data_l, aes(x = Activity_Start_Date, y = Result_Value,
                                            color = Characteristic_Name,
                                            fill = Characteristic_Name))+
    geom_point(pch = 1, color = "black", stroke = 0.2, stat = "identity", show.legend = F) +
    geom_point(data = data_l[data_l$Subject_Taxonomic_Name == "water" , ], 
               aes(x = Activity_Start_Date, y = Result_Value),
               pch = 1, color = "dodgerblue3", stroke = 0.5, size = 1, stat = "identity", show.legend = F) +
    geom_point(data= data_l[data_l$Characteristic_Name == "perfluorooctane sulfonic acid", ], 
               aes(x = Activity_Start_Date, y = Result_Value),
               pch = 1, color = "red3",size = 1, stroke = 0.5, stat = "identity", show.legend = F) +
    theme_bw()+
    ylim(0,250)+
    scale_fill_viridis_d()+
    scale_color_viridis_d()+
    theme(axis.text.x = element_text(angle = 90, face = "italic"))+
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          panel.spacing = unit(0, "lines"))
  
  cowplot::plot_grid(p1.b, p2.b, nrow = 2)
  
  
  # correlation between water and fish tissue for each PFAA
  

# Data summaries ---- 


## Figure 1 ---- 

## Figure 2 ----

