# title: "NRP and JBA new data"
# author: "Abbi Brown (base code for figures) Krista Kraskura (updates)"
# date: jan 6 2025

# *******************************
# Libraries ------ 
library(ggforce)
library(ggfortify)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(tidyverse)
library(here)


# *******************************
# Data ------ 
jba_all<-read.csv(here("./Data/SERDP_report_support/exports/jba_compiled_jan2025.csv"))
wg_all<-read.csv(here("./Data/SERDP_report_support/exports/wg_compiled_jan2025.csv"))

# *******************************
# Colors ------ 
P50<- c("11ClPF3OUdS" = "#FD0022",
 "9ClPF3ONS" = "red4",
 "ADONA" = "black",
 "82FTCA" = "#FB0DE2",
 "82FTS"  = "#00FBE2",
 "82FTUCA" = "#F2C4C7",
 "82diPAP" = "#FE948F", 
 "53FTCA" = "#B60D5D",
 "EtFOSAA"  = "#F2E300",
 "FBSA"  = "grey80",
 "FOSA"   = "#D491FF",
 "42FTS" = "#843573",
 "HFPODA"  = "green",     
"MeFOSA"  = "#AAF597",
"MeFOSAA"  = "#36C2CE",
"NEtFOSA"  = "#0081cf",
"NEtFOSAA" = "yellow",
"NEtFOSE"  = "blue4",
"NFDHA" = "#0D95FC",
"NMeFOSA"  = "#62686A",
"NMeFOSAA" = "grey30",
"73FTCA"  = "#A49949",
"62DiPAP" = "#FB78D1",
"62FTCA"  = "#7A3800",
"62FTS"  = "#BF0DFF",
"62FTUCA" = "#AAF597",
"102FTCA"  = "#843573",
"102FTS"  = "#BEE7FB",
"NMeFOSE" = "#79BE9E",
"PFBA" = "#FF0060",
"PFBS" = "#FFFF80",
"PFDoA"  = "grey19",
"PFDoS"  = "brown",
"PFDS"  = "#005600",
"PFEESA"  = "white",
"PFHpA" = "#00c9a7",
"PFHpS" = "#402E7A",
"PFHxA"  = "#C63C51",
"FHxSA"  = "yellow4",
"PFHxDA" = "#C9DABF",
"PFMBA"  = "#f3c5ff",
"PFMPA"  = "#B5C18E",
"PFNS"  = "#1679A4" ,
"PFODA"  = "#FFFF80",
"diSAmPAP"  = "#E9C874",
 "NMeFOSAA"  = "#88D66C",
"FOSAA"  = "#059212",
 'PFEtCHxS' = '#AA4499',
 'PFHxS' = '#332288', # 
 'PFOS' = "#C4A2D5", # 
 'PFOA' = '#117733',
 'PFNA' = '#88CCEE',
 'PFDA' = 'orange',
"PFOSA" ="#C0D6E8" ,
"PFPeA" = "#E9C874",
"PFPeS"  = '#882255',
"PFPrS" = "#FC1CB1",
"PFTeDA" = "#AFF700",
"PFTrDA" = "#D20062",
"ClPFOS" ="#ECFFE6",
"PFUnA" = "#0D1282",        
"33FTCA" = "#22D9FF", 
 "EtFOSA" = "grey",
 "MeFOSE"= "#FD5C5A",
  "EtFOSE" = "#22D943",
 "PFecHS" = "grey50")




# *******************************
# Figures  -----

# Notes: Martin et al 2003 show that equilibrium of PFAS 
#        uptake happens within 30 days. 

# prep sum PFAS data WG and JBA per sample  ----
# jba (AFRC)
sumJBA<-rbind(jba_all) %>% 
  filter(Media == "water") %>% 
  mutate(SampleID = factor(SampleID)) %>% 
  dplyr::group_by(Media, SampleDate, Spp_Loc, Seasonality) %>% 
  dplyr::summarise(SumPFAS = sum(value, na.rm = T)) # ng/L
sumJBAsed<-rbind(jba_all) %>% 
  filter(Media == "sediment") %>% 
  mutate(SampleID = factor(SampleID)) %>% 
  dplyr::group_by(Media, SampleDate, Spp_Loc, Seasonality) %>% 
  dplyr::summarise(SumPFAS = sum(value, na.rm = T)) # ng/L
sumJBAmuscle<-rbind(jba_all) %>% 
  filter(Media == "muscle") %>% 
  mutate(SampleID = factor(SampleID)) %>% 
  dplyr::group_by(Media, SampleDate, SampleID, Spp_Loc, Seasonality) %>% 
  dplyr::summarise(SumPFAS = sum(value, na.rm = T)) # ng/kg
sumJBAliver<-rbind(jba_all) %>% 
  filter(Media == "liver") %>% 
  mutate(SampleID = factor(SampleID)) %>% 
  dplyr::group_by(Media, SampleDate, SampleID, Spp_Loc,Seasonality) %>% 
  dplyr::summarise(SumPFAS = sum(value, na.rm = T)) # ng/kg
sumJBAfish<-rbind(jba_all) %>% 
  filter(Media == "whole body") %>% 
  mutate(SampleID = factor(SampleID)) %>% 
  dplyr::group_by(Media, SampleDate, SampleID, Spp_Loc, Seasonality) %>% 
  dplyr::summarise(SumPFAS = sum(value, na.rm = T)) # ng/kg

sumJBA_all<-rbind(sumJBA, 
      sumJBAfish[, -3],
      sumJBAliver[, -3], 
      sumJBAmuscle[, -3])

sumJBA_fish<-rbind(
      sumJBAfish[, -3],
      sumJBAliver[, -3], 
      sumJBAmuscle[, -3])

# wg (NRP) 
sumWG<-rbind(wg_all) %>% 
  filter(Media == "water") %>% 
  mutate(SampleID = factor(SampleID)) %>% 
  dplyr::group_by(Media, SampleDate, Spp_Loc, Seasonality) %>% 
  dplyr::summarise(SumPFAS = sum(value, na.rm = T)) # ng/L
sumWGmuscle<-rbind(wg_all) %>% 
  filter(Media == "muscle") %>% 
  mutate(SampleID = factor(SampleID)) %>% 
  dplyr::group_by(Media, SampleDate, SampleID, Spp_Loc, Seasonality) %>% 
  dplyr::summarise(SumPFAS = sum(value, na.rm = T)) # ng/kg
sumWGliver<-rbind(wg_all) %>% 
  filter(Media == "liver") %>% 
  mutate(SampleID = factor(SampleID)) %>% 
  dplyr::group_by(Media, SampleDate, SampleID, Spp_Loc, Seasonality) %>% 
  dplyr::summarise(SumPFAS = sum(value, na.rm = T)) # ng/kg
sumWGfish<-rbind(wg_all) %>% 
  filter(Media == "whole body") %>% 
  mutate(SampleID = factor(SampleID)) %>% 
  dplyr::group_by(Media, SampleDate, SampleID, Spp_Loc, Seasonality) %>% 
  dplyr::summarise(SumPFAS = sum(value, na.rm = T)) # ng/kg

sumWG_all<-rbind(sumWG, 
      sumWGfish[, -3],
      sumWGliver[, -3], 
      sumWGmuscle[, -3])

sumWG_fish<-rbind( 
      sumWGfish[, -3],
      sumWGliver[, -3], 
      sumWGmuscle[, -3])

sumWG_fish$System <- "NRP"
sumJBA_fish$System <- "AFBC"
sumFish<-rbind(sumJBA_fish, sumWG_fish)


# *******************************
# timeline plot ----
p1<-ggplot(sumWG, aes(x=as.Date(SampleDate), y=SumPFAS,
                  color = Seasonality)) +
    geom_point(data = sumJBA,
               aes(x=as.Date(SampleDate),
                   y=SumPFAS,
                   color = Seasonality), pch = 21, size = 1, stroke = 0.3) +
    geom_point(size=1, pch = 19)+
    ylab("Sum PFAS Conc. (ng/L)") +
    xlab("")+
    scale_color_manual(values = c("blue4", "green4", "red4"))+
    scale_fill_manual(values = c("blue4", "green4", "red4"))+
    scale_x_date(date_labels = "%b %d", breaks = "1 month", 
                limits = c(as.Date("2020-10-01"), as.Date("2021-09-01")))+
    theme_bw()+
    ylim(0, 10000)+
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          panel.grid = element_blank(),
          legend.position = "none")+
  ggtitle("NRP - closed, AFBC - open")
p1



# *******************************
# WATER  ------
# *******************************
# *******************************
# timeline plot of different PFAS 
p1wg<-wg_all %>% 
  filter(Media == "water") %>% 
  group_by(SampleDate, Analyte, Seasonality) %>% 
  summarize(mean_val_day = mean(value, na.rm = TRUE)) %>%  
  ggplot(aes(x=as.Date(SampleDate),
               y= mean_val_day,
               fill = Analyte,
               color = Analyte, 
               group = interaction(Analyte))) +
    geom_hline(yintercept = 1000, color = "grey", linetype = "dashed", 
               linewidth = 0.1)+
    geom_bar(position='stack', width = 1.5, stat='identity') +
    ylab("PFAS in water (ng/L)") +
    xlab("")+
    coord_flip()+
    # geom_smooth(size = 0.5)+
    # geom_smooth(data = sumJBA, size = 0.5)+
    scale_color_manual(values = P50)+
    scale_fill_manual(values = P50)+
    scale_x_date(date_labels = "%b %d %Y", breaks = "1 month", 
                limits = c(as.Date("2020-10-05"), as.Date("2021-09-01")))+
    theme_bw()+
    ylim(0, 8000)+
    theme(axis.text.x = element_text(angle = 90, hjust=1),
          panel.grid = element_blank(),
          legend.position = "none")
  # ggtitle("NRP")
p1wg

ggsave(p1wg, filename = "./Figures/SERDP/comparison_barStackNRP.png",
       height = 5, width = 2.5, units = "in")


# *******************************
# *******************************
# *******************************
# timeline plot of different PFAS
p1jba<-jba_all %>% 
  filter(Media == "water") %>% 
  group_by(SampleDate, Analyte, Seasonality) %>% 
  summarize(mean_val_day = mean(value, na.rm = TRUE)) %>%  
  ggplot(aes(x=as.Date(SampleDate),
               y= mean_val_day,
               fill = Analyte,
               color = Analyte, 
               group = interaction(Analyte))) +
    geom_hline(yintercept = 1000, color = "grey", linetype = "dashed", 
               linewidth = 0.1)+
    geom_bar(position='stack', width = 1.5, stat='identity') +
    # facet_grid(.~Seasonality)+
    ylab("PFAS in water (ng/L)") +
    xlab("")+
    coord_flip()+
    # geom_smooth(size = 0.5)+
    # geom_smooth(data = sumJBA, size = 0.5)+
    scale_color_manual(values = P50)+
    scale_fill_manual(values = P50)+
    scale_x_date(date_labels = "%b %d %Y",
                 breaks = "1 month", 
                 limits = c(as.Date("2020-10-05"), as.Date("2021-09-01")))+
    theme_bw()+
    ylim(0, 8000)+
    theme(axis.text.x = element_text(angle = 90, hjust=1),
          panel.grid = element_blank(),
          legend.position = "none")
  # ggtitle("AFBC") 
ggsave(p1jba, filename = "./Figures/SERDP/comparison_barStackAFBC.png",
       height = 5, width = 2.5, units = "in")



# *******************************
# *******************************
# *******************************
# FISH PERCENT PFOS
wg_pct_fish<-wg_all %>% 
  filter(Media == "whole body") %>% 
  group_by(SampleID, SampleDate) %>%
  mutate(sumPFAS_sample = sum(value, na.rm = TRUE)) %>%
  group_by(SampleID, Analyte, SampleDate) %>%
  mutate(pctVal = value/sumPFAS_sample * 100) %>% 
  select(Analyte, SampleID, Spp_Loc,
         SampleDate, Media, System,
         Seasonality, value, pctVal, sumPFAS_sample)%>% 
  group_by(SampleDate, Analyte) %>% 
  summarize(meanPctVal = mean(pctVal),
            sd = sd(pctVal)) %>% 
  filter(Analyte == "PFOS")

jba_pct_fish<-jba_all %>% 
  filter(Media == "whole body") %>% 
  group_by(SampleID, SampleDate) %>%
  mutate(sumPFAS_sample = sum(value, na.rm = TRUE)) %>%
  group_by(SampleID, Analyte, SampleDate) %>%
  mutate(pctVal = value/sumPFAS_sample * 100) %>% 
  select(Analyte, SampleID, Spp_Loc,
         SampleDate, Media, System,
         Seasonality, value, pctVal, sumPFAS_sample) %>% 
  group_by(SampleDate, Analyte) %>% 
  summarize(meanPctVal = mean(pctVal),
            sd = sd(pctVal)) %>% 
  filter(Analyte == "PFOS") 

# *******************************
# WATER PERCENT PFOS
wg_pct_water<-wg_all %>% 
  filter(Media == "water") %>% 
  group_by(SampleID, SampleDate) %>%
  mutate(sumPFAS_sample = sum(value, na.rm = TRUE)) %>%
  group_by(SampleID, Analyte, SampleDate) %>%
  mutate(pctVal = value/sumPFAS_sample * 100) %>% 
  select(Analyte, SampleID, Spp_Loc,
         SampleDate, Media, System,
         Seasonality, value, pctVal, sumPFAS_sample)

jba_pct_water<-jba_all %>% 
  filter(Media == "water") %>% 
  group_by(SampleID, SampleDate) %>%
  mutate(sumPFAS_sample = sum(value, na.rm = TRUE)) %>%
  group_by(SampleID, Analyte, SampleDate) %>%
  mutate(pctVal = value/sumPFAS_sample * 100) %>% 
  select(Analyte, SampleID, Spp_Loc,
         SampleDate, Media, System,
         Seasonality, value, pctVal, sumPFAS_sample)

p2PFOS_jba<-jba_pct_water %>% 
  group_by(SampleDate, Analyte) %>% 
  summarize(meanPctVal = mean(pctVal),
            sd = sd(pctVal)) %>% 
  filter(Analyte == "PFOS") %>% 
  ggplot(aes(x = as.Date(SampleDate), y = meanPctVal,
             color = Analyte,
             fill = Analyte))+
  geom_hline(yintercept = c(25, 50,75, 95), color = "grey30", linetype = "dashed", 
               linewidth = 0.1)+
  geom_errorbar(aes(ymin = meanPctVal - sd, ymax = meanPctVal + sd),
                color = "black", size = 0.3, width = 0.1)+
  geom_point(pch = 21, color= "black", stroke = 0.2, size = 2)+
  geom_errorbar(data = jba_pct_fish, 
                aes(ymin = meanPctVal - sd, ymax = meanPctVal + sd),
                color = "black", size = 0.3, width = 0.1)+
  geom_point(data = jba_pct_fish, pch = 21, color= "black", size = 3, stroke = 0.7)+  
  scale_color_manual(values = P50)+
  scale_fill_manual(values = P50)+
  scale_x_date(date_labels = "%b %d %Y",
               breaks = "1 month", 
               limits = c(as.Date("2020-10-05"), as.Date("2021-09-01")))+
  coord_flip()+
  theme_bw()+
  ylim(0, 100)+
  ylab("Percent")+
  theme(axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")
p2PFOS_jba
ggsave(p2PFOS_jba, filename = "./Figures/SERDP/comparison_PFOS_AFBC.png",
       height = 5, width = 4.5, units = "in")


p2PFOS_wg<-wg_pct_water %>% 
  group_by(SampleDate, Analyte) %>% 
  summarize(meanPctVal = mean(pctVal),
            sd = sd(pctVal)) %>% 
  filter(Analyte == "PFOS") %>% 
  ggplot(aes(x = as.Date(SampleDate), y = meanPctVal,
             color = Analyte,
             fill = Analyte))+
  geom_hline(yintercept = c(25, 50,75, 95), color = "grey30", linetype = "dashed", 
               linewidth = 0.1)+
  geom_errorbar(aes(ymin = meanPctVal - sd, ymax = meanPctVal + sd),
                color = "black", size = 0.3, width = 0.1)+
  geom_point(pch = 21, color= "black", stroke = 0.2, size = 2)+
  geom_errorbar(data = wg_pct_fish, 
                aes(ymin = meanPctVal - sd, ymax = meanPctVal + sd),
                color = "black", size = 0.3, width = 0.1)+
  geom_point(data = wg_pct_fish, pch = 21, color= "black", size = 3, stroke = 0.7)+  
  scale_color_manual(values = P50)+
  scale_fill_manual(values = P50)+
  scale_x_date(date_labels = "%b %d %Y",
               breaks = "1 month", 
               limits = c(as.Date("2020-10-05"), as.Date("2021-09-01")))+
  coord_flip()+
  theme_bw()+
  ylim(0, 100)+
  ylab("Percent")+
  theme(axis.title.y = element_blank(),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")
p2PFOS_wg
ggsave(p2PFOS_wg, filename = "./Figures/SERDP/comparison_PFOS_NRP.png",
       height = 5, width = 4.5, units = "in")


# Percent of the major PFAS > 5 % in each system 
wg.pct<-wg_pct_water %>% 
  group_by(SampleDate, Analyte) %>% 
  summarize(meanPctVal = mean(pctVal),
            sd = sd(pctVal)) %>% 
  filter(meanPctVal>5) %>% 
  ggplot(aes(x = as.Date(SampleDate), y = meanPctVal,
             color = Analyte,
             fill = Analyte))+
  # geom_hline(yintercept = c(25, 50,75, 95), color = "grey30", linetype = "dashed", 
               # linewidth = 0.1)+
  geom_errorbar(aes(ymin = meanPctVal - sd, ymax = meanPctVal + sd),
                size = 0.3, width = 0.1)+
  geom_point(pch = 21, stroke = 0.2, size = 3, alpha = 0.8)+
  scale_color_manual(values = P50)+
  scale_fill_manual(values = P50)+
  scale_x_date(date_labels = "%b %d",
               breaks = "1 month", 
               limits = c(as.Date("2020-10-05"), as.Date("2021-09-01")))+
  ylim(0,100)+
  theme_pubr()+
  theme(legend.position = "top", 
        legend.title = element_blank())+
  labs(y = "Percent of SumPFAS", x = "")+
  coord_flip()

jba.pct<-jba_pct_water %>% 
  group_by(SampleDate, Analyte) %>% 
  summarize(meanPctVal = mean(pctVal),
            sd = sd(pctVal)) %>% 
  filter(meanPctVal>5) %>% 
  ggplot(aes(x = as.Date(SampleDate), y = meanPctVal,
             color = Analyte,
             fill = Analyte))+
  geom_errorbar(aes(ymin = meanPctVal - sd, ymax = meanPctVal + sd),
                size = 0.3, width = 0.1)+
  geom_point(pch = 21, stroke = 0.2, size = 3, alpha = 0.8)+
  scale_color_manual(values = P50)+
  scale_fill_manual(values = P50)+
  scale_x_date(date_labels = "%b %d",
               breaks = "1 month", 
               limits = c(as.Date("2020-10-05"), as.Date("2021-09-01")))+
  theme_pubr()+
  ylim(0,100)+
  theme(legend.position = "top", 
        legend.title = element_blank())+
  labs(y = "Percent of SumPFAS", x = "")+
  coord_flip()

ggsave(cowplot::plot_grid( wg.pct, jba.pct,
                          labels = c("NRP", "AFBC"), ncol =2),
       filename = "./Figures/SERDP/above5percent_sumPFAS.png",
       height = 5, width = 10, units = "in")



# *******************************
# FISH BARPLOTS ----------
# *******************************
# *******************************
# season plot of different PFAS  in fish muscle
p1wg.f<-wg_all %>% 
  filter(Media == "whole body") %>% 
  group_by(SampleDate, Spp_Loc, Analyte, Seasonality) %>% 
  summarize(mean_val_day = mean(value, na.rm = TRUE), n = n()) %>%  
  mutate(spp_season = paste(Spp_Loc, Seasonality)) %>% 
  ggplot(aes(x=Spp_Loc,
               y= mean_val_day,
               fill = Analyte,
               color = Analyte)) +
    geom_hline(yintercept = 1000, color = "grey", linetype = "dashed", 
               linewidth = 0.1)+
    geom_bar(position='stack', width = 0.8, stat='identity') +
    ylab("PFAS in fish muscle (ng/g)") +
    facet_wrap(.~Seasonality, ncol = 1)+
    coord_flip()+
    xlab("")+
    scale_color_manual(values = P50)+
    scale_fill_manual(values = P50)+
    theme_bw()+
    ylim(0, 2100)+
    theme(axis.text.x = element_text(angle = 90, hjust=1),
          panel.grid = element_blank(),
          legend.position = "none")+
  ggtitle("NRP")
p1wg.f

ggsave(p1wg.f, filename = "./Figures/SERDP/comparison_barStackNRP_fish.png",
       height = 2, width = 2.4, units = "in")

# timeline plot of different PFAS 
p1jba.f<-jba_all %>% 
  filter(Media == "whole body") %>% 
  group_by(SampleDate, Spp_Loc, Analyte, Seasonality) %>% 
  summarize(mean_val_day = mean(value, na.rm = TRUE), n = n()) %>%  
  mutate(spp_season = paste(Spp_Loc, Seasonality)) %>% 
  ggplot(aes(x=Spp_Loc,
               y= mean_val_day,
               fill = Analyte,
               color = Analyte)) +
    geom_hline(yintercept = 1000, color = "grey", linetype = "dashed",
               linewidth = 0.1)+
    geom_bar(position='stack', width = 0.8, stat='identity') +
    ylab("PFAS in fish tissue (ng/g)") +
    facet_wrap(.~Seasonality, ncol = 1)+
    coord_flip()+
    xlab("")+
    scale_color_manual(values = P50)+
    scale_fill_manual(values = P50)+
    theme_bw()+
    ylim(0, 2100)+
    theme(axis.text.x = element_text(angle = 90, hjust=1),
          panel.grid = element_blank(),
          legend.position = "none")+
  ggtitle("AFBC")
p1jba.f

ggsave(p1jba.f, filename = "./Figures/SERDP/comparison_barStackAFBC_fish.png",
       height = 5, width = 2.5, units = "in")


### BAF -----
# take any fish and any day in the season 
# what is BAF error and predictability?

wg_all %>% 
  select(Analyte, SampleDate, SampleID,
         Spp_Loc, Weight, value,
         Media, Seasonality, Chain_len,
         Acronym2, System) 

BAFcalc<-function(data, Analyte, Season, tissue, niter){
  
  f_sum<-data[c(data$Analyte == Analyte &  
         data$Media == tissue &  
         data$Seasonality == Season),]
  
  w_sum<-data[c(data$Analyte == Analyte &  
         data$Media == "water" &  
         data$Seasonality == Season),]


  for(i in 1:niter){
    # random sample row 
    f_r<-sample(nrow(f_sum), 1)
    w_r<-sample(nrow(w_sum), 1)
    new_row<-f_sum[f_r,]
    new_row$water_value_ngML<-w_sum[w_r, "value"]/1000
    new_row$logBAF<-NA
  
    if (i ==1){
      new_row$logBAF<-log10(new_row$value/new_row$water_value_ngML)
      d<-new_row
    }else{
      new_row$logBAF<-log10(new_row$value/new_row$water_value_ngML)
      d<-rbind(d, new_row)
    }
    if(i == niter){
      for(j in 1:nrow(d)){
        d$logBAFquant[j]<-ecdf(c(d$logBAF))(d$logBAF[j])
        if(d$logBAFquant[j] > 0.5) {
          d$logBAFquant[j] <- 1-d$logBAFquant[j]
        }
      }
    }
  }
  
  return(d)
}


d_wg_summer<-BAFcalc(data = wg_all,
        Analyte = "PFOS",
        Season = "summer",
        tissue = "whole body",
        niter = 10000)

d_wg_winter<-BAFcalc(data = wg_all,
        Analyte = "PFOS",
        Season = "fall/winter",
        tissue = "whole body",
        niter = 10000)

d_jba_summer<-BAFcalc(data = jba_all,
        Analyte = "PFOS",
        Season = "summer",
        tissue = "whole body",
        niter = 10000)

d_jba_spring<-BAFcalc(data = jba_all,
        Analyte = "PFOS",
        Season = "spring",
        tissue = "whole body",
        niter = 10000)

ggplot(data = d_wg_summer, aes(x = logBAF, y = logBAFquant*100))+
  geom_point(pch = 21, size = 1)+
  geom_point(data = d_wg_winter,mapping = aes(x = logBAF, y = logBAFquant*100),
             pch = 21, size = 1, color = "grey")+
  geom_point(data = d_jba_summer,mapping = aes(x = logBAF, y = logBAFquant*100),
             pch = 23, size = 1, color = "red3")+
  geom_point(data = d_jba_spring,mapping = aes(x = logBAF, y = logBAFquant*100),
             pch = 23, size = 1, color = "pink")+
  theme_pubr()+
  xlim(0,5.2)+
  facet_wrap(.~System, nrow = 2)+
  ylab("Percentile")+
  xlab(expression(log[10]~BAF~(ng/g)~w*w))


# same with PFHxS
d_wg_summer<-BAFcalc(data = wg_all,
        Analyte = "PFHxS",
        Season = "summer",
        tissue = "whole body",
        niter = 10000)

d_wg_winter<-BAFcalc(data = wg_all,
        Analyte = "PFHxS",
        Season = "fall/winter",
        tissue = "whole body",
        niter = 10000)

d_jba_summer<-BAFcalc(data = jba_all,
        Analyte = "PFHxS",
        Season = "summer",
        tissue = "whole body",
        niter = 10000)

d_jba_spring<-BAFcalc(data = jba_all,
        Analyte = "PFHxS",
        Season = "spring",
        tissue = "whole body",
        niter = 10000)

ggplot(data = d_wg_summer, aes(x = logBAF, y = logBAFquant*100))+
  geom_point(pch = 21, size = 1)+
  geom_point(data = d_wg_winter,mapping = aes(x = logBAF, y = logBAFquant*100),
             pch = 21, size = 1, color = "grey")+
  geom_point(data = d_jba_summer,mapping = aes(x = logBAF, y = logBAFquant*100),
             pch = 23, size = 1, color = "red3")+
  geom_point(data = d_jba_spring,mapping = aes(x = logBAF, y = logBAFquant*100),
             pch = 23, size = 1, color = "pink")+
  theme_pubr()+
  xlim(0,5.2)+
  facet_wrap(.~System, nrow = 2)+
  ylab("Percentile")+
  xlab(expression(log[10]~BAF~(ng/g)~w*w))


# *******************************************
# *******************************************
# *******************************************
# sample dates
 # WG/NRP   as.Date("2020-11-10")
            # as.Date("2021-04-28")
            # as.Date("2021-08-19")

# JBA/AFBC as.Date("2021-04-28")
         # as.Date("2021-07-16")

# WATER ----
d_wg_wat_baf <- d_wg_wat %>%
  dplyr::select(SampleID, SampleDate, all_of(wgPFAS)) %>% 
  pivot_longer(cols = all_of(wgPFAS),
         names_to = "Analyte") %>%
  mutate(Media = "water",
         SampleDate = strptime(SampleDate, format = "%m/%d/%Y %H:%M")) %>% 
  mutate(Seasonality = if_else(SampleDate > as.Date("2021-06-01"), "summer", 
                              if_else(SampleDate < as.Date("2020-11-20"), "fall/winter",
                                      "spring"))) %>% 
  filter(SampleDate == "2020-11-10" |
         SampleDate == "2021-04-28" |
         SampleDate == "2021-08-19") %>% 
  group_by(Seasonality, Analyte) %>% 
  summarise(meanPFASw = mean(value, na.rm = TRUE))%>% 
  mutate(Seas_Analyte = paste(Seasonality, Analyte, sep="-"))

d_jba_wat_baf <- d_jba_wat %>%
  dplyr::select(SampleID, SampleDate, all_of(jbaPFAS)) %>% 
  pivot_longer(cols = all_of(jbaPFAS),
         names_to = "Analyte") %>%
  mutate(Media = "water",
         SampleDate = strptime(SampleDate, format = "%m/%d/%Y %H:%M")) %>% 
  mutate(Seasonality = if_else(SampleDate > as.Date("2021-06-10"), "summer", 
                              if_else(SampleDate < as.Date("2021-01-25"), "fall/winter", "spring"))) %>% 
  filter(SampleDate == "2021-04-28" |
         SampleDate == "2021-07-16") %>% 
  group_by(Seasonality, Analyte) %>% 
  summarise(meanPFASw = mean(value, na.rm = TRUE))%>% 
  mutate(Seas_Analyte = paste(Seasonality, Analyte, sep="-"))

d_wg_wat_baf<-merge(d_analytes, d_wg_wat_baf)
d_jba_wat_baf<-merge(d_analytes, d_jba_wat_baf)

# FISH TISSUES tissues ------
d_wg_biota_baf <- d_wg_biota %>%
  dplyr::select(Species, SampleID, SampleDate, Tissue, 
                all_of(wgPFAS)) %>% 
  pivot_longer(cols = all_of(wgPFAS),
         names_to = "Analyte") %>%
  mutate(Media = Tissue,
         SampleDate = strptime(SampleDate, format = "%m/%d/%Y %H:%M")) %>% 
  mutate(Seasonality = if_else(SampleDate > as.Date("2021-06-01"), "summer", 
                              if_else(SampleDate < as.Date("2020-11-20"), "fall/winter", "spring"))) %>% 
  group_by(Seasonality, Analyte, Species, Media) %>% 
  summarise(meanPFASbiota = mean(value, na.rm = TRUE))%>% 
  mutate(Seas_Analyte = paste(Seasonality, Analyte, sep="-"), 
         Media = tolower(Media)) 
# %>% 
  # filter(Media == "muscle" | Media == "whole body")

d_jba_biota_baf <- d_jba_biota %>%
  dplyr::select(CommonName, SampleID, SampleDate, Tissue,
               all_of(jbaPFAS)) %>% 
  pivot_longer(cols = all_of(jbaPFAS),
         names_to = "Analyte") %>%
  mutate(Media = Tissue,
         SampleDate = strptime(SampleDate, format = "%m/%d/%Y %H:%M")) %>% 
  mutate(Seasonality = if_else(SampleDate > as.Date("2021-06-10"), "summer", 
                              if_else(SampleDate < as.Date("2021-01-25"), "fall/winter", "spring"))) %>% 
  group_by(Seasonality, Analyte, CommonName, Media) %>% 
  summarise(meanPFASbiota = mean(value, na.rm = TRUE)) %>% 
  mutate(Seas_Analyte = paste(Seasonality, Analyte, sep="-"), 
         Media = tolower(Media)) 
# %>% 
  # filter(Media == "muscle" | Media == "whole body")

d_jba_baf<-merge(d_jba_biota_baf, d_jba_wat_baf[, c("meanPFASw", "Seas_Analyte",  "Chain_len", "Acronym2")], by = "Seas_Analyte")
d_wg_baf<-merge(d_wg_biota_baf, d_wg_wat_baf[, c("meanPFASw", "Seas_Analyte", "Chain_len", "Acronym2")], by = "Seas_Analyte")


# Figures SAME DAY BAFS ONLY (different than the ones in CompleteDatasets_rerun file-------
jba_baf<-d_jba_baf %>% 
  mutate(logBAF = log10(meanPFASbiota/ c(meanPFASw/1000)), 
         Chain_len = as.numeric(Chain_len)) %>% 
  mutate(Analyte = fct_reorder(Analyte, Chain_len)) %>%
  ggplot(aes(x = Analyte, y = logBAF, 
           group = interaction(CommonName, Media),
           color = CommonName,
           fill = CommonName,
           shape = Media)) +
  # geom_boxplot(alpha = 0.5)+
  geom_point(size=2,  alpha = 1)+
  scale_shape_manual(values = c(21, 22, 23, 24))+
  scale_color_futurama()+
  scale_fill_futurama()+
  facet_grid(Seasonality~Acronym2, scales = "free")+
  xlab("")+
  ylab(expression(log[10]~BAF~(ng/g)~w*w))+ 
  ggtitle("AFBC")+
  geom_hline(yintercept = 4, 
             color= "grey",
             linewidth = 0.2,
             linetype = "dashed")+
  # geom_smooth(method = "lm", se = FALSE)+
  theme_bw()+
  guides(color=guide_legend(ncol=1))+
  theme(legend.direction = "vertical", 
        legend.box = "vertical",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c("right"),
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust =1))
ggsave(jba_baf,filename = "./Figures/SERDP/JBA_BAF_same_dayonly.png",
       width = 13, height = 6.5, units = "in")


wg_baf<-d_wg_baf %>% 
  mutate(logBAF = log10(meanPFASbiota/ c(meanPFASw/1000)), 
         Chain_len = as.numeric(Chain_len)) %>% 
  mutate(Analyte = fct_reorder(Analyte, Chain_len)) %>%
  ggplot(aes(x = Analyte, y = logBAF, 
           group = interaction(Species, Media),
           color = Species,
           fill = Species,
           shape = Media, 
           label = Chain_len)) +
  geom_point(size=2, alpha = 1)+
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  scale_color_frontiers()+
  scale_fill_frontiers()+
  facet_grid(Seasonality~Acronym2, scales = "free")+
  xlab("")+
  ylab(expression(log[10]~BAF~(ng/g)~w*w))+ 
  ggtitle("NRP")+
  geom_hline(yintercept = 4, 
             color= "grey",
             linewidth = 0.2,
             linetype = "dashed")+
  # geom_smooth(method = "lm", se = FALSE)+
  theme_bw()+
  # geom_text(mapping = aes(y = -1, x = Analyte), color = "black")+
  # ylim(-1.5, 4.5)+
  guides(color=guide_legend(ncol=1))+
  theme(legend.direction = "vertical", 
        legend.box = "vertical",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = "right",
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust =1))
ggsave(wg_baf,filename = "./Figures/SERDP/WG_BAF_same_dayonly.png",
       width = 13, height = 6.5, units = "in")
  


# ************************
# ************************
# ************************
# The higher the % in the water the higher BAF? 
d_jba_baf0<-d_jba_baf %>% 
  mutate(logBAF = log10(meanPFASbiota/ c(meanPFASw/1000)), 
         Chain_len = as.numeric(Chain_len)) %>% 
  mutate(Analyte = fct_reorder(Analyte, Chain_len))

d_jba_baf_pct_day<-d_jba_baf0 %>% 
  filter (! duplicated(Seas_Analyte)) %>% 
  select(Seas_Analyte, Seasonality, Analyte, meanPFASw) %>% 
  group_by(Seasonality) %>% 
  mutate(sumPFAS = sum(meanPFASw, na.rm = TRUE), 
         pctVal = meanPFASw/sumPFAS * 100)
  
d_wg_baf0<-d_wg_baf %>% 
  mutate(logBAF = log10(meanPFASbiota/ c(meanPFASw/1000)), 
         Chain_len = as.numeric(Chain_len)) %>% 
  mutate(Analyte = fct_reorder(Analyte, Chain_len))

d_wg_baf_pct_day<-d_wg_baf0 %>% 
  filter (! duplicated(Seas_Analyte)) %>% 
  select(Seas_Analyte, Seasonality, Analyte, meanPFASw) %>% 
  group_by(Seasonality) %>% 
  mutate(sumPFAS = sum(meanPFASw, na.rm = TRUE), 
         pctVal = meanPFASw/sumPFAS * 100)

d_jba_baf1<-merge(d_jba_baf0, d_jba_baf_pct_day, by = "Seas_Analyte")
d_wg_baf1<-merge(d_wg_baf0, d_wg_baf_pct_day, by = "Seas_Analyte")

# 
ggplot(d_wg_baf1,
       aes(pctVal, logBAF,color = Analyte.y, 
                      ))+
  geom_point(pch = 19, size=0.6)+
  geom_point(d_jba_baf1,
             mapping = aes(pctVal, logBAF,color = Analyte.y),
             shape = 21, size =0.6)+
  scale_color_manual(values = P50)+
  scale_fill_manual(values = P50)+
  theme_pubr()+
  facet_wrap(~Acronym2, nrow = 3)+
  # theme(legend.position = "none")+
  xlab("percent PFOS of total PFAS in water")
  
  
  
d_wg_baf1 %>%
  filter(Analyte.x == "PFOS" | Analyte.x == "PFHxS") %>% 
  pivot_wider(names_from = station, values_from = seen)



