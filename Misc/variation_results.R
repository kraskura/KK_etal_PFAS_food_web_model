AO<-read.csv("./Data/outputs/mcmc_replicate_runs_observedDiet 2024-01-25 .csv")
library(ggalt)
library(ggrepel)

p_lines <- data.frame('x' = c(0, 7),
                      'y' = c(0, 7))
p_lines$y1 <- p_lines$x-1
p_lines$y2 <- p_lines$x+1
p_lines$y3 <- p_lines$x-log10(2)
p_lines$y4 <- p_lines$x+log10(2)

PFAA_cols <- c('PFHxS' = '#332288', # 
             'PFOS' = '#882255', # 
             'PFOA' = '#117733',
             'PFNA' = '#88CCEE',
             'PFDA' = 'orange',
             'PFUA' = '#AA4499')
# set percent error: 
AO<-AO[!c(AO$PFAA == "PFUA" | AO$PFAA == "PFDA"),]
AO$system<-paste(AO$loc, AO$Seasonality, sep =" ")

AO$spp_PFAA<-paste(AO$SppAlias, AO$PFAA, sep =" ")

ggplot(data = AO, aes(x = log_ngkg_re, y = Obs_logngkg, color = temperature))+
    facet_grid(.~system)+
    geom_point(pch = 19, size = 0.01)+
    theme_bw()+
    ggtitle(label = 'temperature')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')#+

# mass 
ggplot(data = AO[AO$system == "JBA Spring", ], aes(x = log_ngkg_re, y = Obs_logngkg, color = WB))+
    facet_grid(.~system)+
    geom_point(pch = 19, size = 0.01)+
    theme_bw()+
    ggtitle(label = 'Body mass')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')#+

ggplot(data = AO[AO$system == "JBA Summer", ], aes(x = log_ngkg_re, y = Obs_logngkg, color = WB))+
    facet_grid(.~system)+
    geom_point(pch = 19, size = 0.01)+
    theme_bw()+
    ggtitle(label = 'Body mass')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')#+

ggplot(data = AO[AO$system == "WG Fall", ], aes(x = log_ngkg_re, y = Obs_logngkg, color = WB))+
    facet_grid(.~system)+
    geom_point(pch = 19, size = 0.01)+
    theme_bw()+
    ggtitle(label = 'Body mass')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')#+

ggplot(data = AO, aes(x = log_ngkg_re, y = Obs_logngkg, color = env_DO))+
    facet_grid(.~system)+
    geom_point(pch = 19, size = 0.01)+
    theme_bw()+
    ggtitle(label = 'Dissolved O2')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')#+

ggplot(data = AO, aes(x = log_ngkg_re, y = Obs_logngkg, color = C_s))+
    facet_grid(.~system)+
    geom_point(pch = 19, size = 0.01)+
    theme_bw()+
    ggtitle(label = 'PFAA in sediment')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')#+

ggplot(data = AO, aes(x = log_ngkg_re, y = Obs_logngkg, color = C_WTO))+
    facet_grid(.~system)+
    geom_point(pch = 19, size = 0.01)+
    theme_bw()+
    ggtitle(label = 'PFAA in water')+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')#+


ggplot(data = AO, aes(x = log_ngkg_re, y = Obs_logngkg,
                      color = spp_PFAA))+
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
    facet_grid(.~system)+
    # geom_point(pch = 19, size = 0.01)+
    theme_bw()+
    geom_encircle(aes(fill = spp_PFAA, color = spp_PFAA),
                  s_shape = 1, expand = 0,
                  alpha = 0.3,  show.legend = F)+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')#+


# Bgl ------
# Species by temperature
ggplot(data = AO[AO$SppAlias == "Bgl",], aes(x = log_ngkg_re, y = Obs_logngkg,
                      color = temperature))+
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
    facet_grid(.~system)+
    scale_color_viridis_c()+
    theme_classic()+
    geom_encircle(aes(fill = SppAlias, color = SppAlias),
                  s_shape = 1, expand = 0, color = "black", fill = "grey",
                  alpha = 0.3,  show.legend = F)+
    geom_point(pch = 19, size = 0.2)+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    ggtitle("Bgl")#+

# Fal ------
ggplot(data = AO[AO$SppAlias == "Fal",], aes(x = log_ngkg_re, y = Obs_logngkg,
                      color = temperature))+
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
    facet_grid(.~system)+
    geom_encircle(aes(fill = SppAlias, color = SppAlias),
                  s_shape = 1, expand = 0, color = "black", fill = "grey",
                  alpha = 0.3,  show.legend = F)+
    geom_point(pch = 19, size = 0.2)+
    scale_color_viridis_c()+
    theme_classic()+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    ggtitle("Fal")#+

# Pum ------
ggplot(data = AO[AO$SppAlias == "Pum",], aes(x = log_ngkg_re, y = Obs_logngkg,
                      color = temperature))+
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
    facet_grid(.~system)+
    scale_color_viridis_c()+
    theme_classic()+
    geom_encircle(aes(fill = SppAlias, color = SppAlias),
                  s_shape = 1, expand = 0, color = "black", fill = "grey",
                  alpha = 0.3,  show.legend = F)+
    geom_point(pch = 19, size = 0.2)+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    ggtitle("Pum")#+

# Bas ------
ggplot(data = AO[AO$SppAlias == "Bas",], aes(x = log_ngkg_re, y = Obs_logngkg,
                      color = temperature))+
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
    facet_grid(.~system)+
    geom_encircle(aes(fill = SppAlias, color = SppAlias),
                  s_shape = 1, expand = 0, color = "black", fill = "grey",
                  alpha = 0.3,  show.legend = F)+
    geom_point(pch = 19, size = 0.2)+
    scale_color_viridis_c()+
    theme_classic()+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    ggtitle("Bas")#+


# Swa ------
ggplot(data = AO[AO$SppAlias == "Swa",], aes(x = log_ngkg_re, y = Obs_logngkg,
                      color = temperature))+
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
    facet_grid(.~system)+
    scale_color_viridis_c()+
    theme_classic()+
    geom_encircle(aes(fill = SppAlias, color = SppAlias),
                  s_shape = 1, expand = 0, color = "black", fill = "grey",
                  alpha = 0.3,  show.legend = F)+
    geom_point(pch = 19, size = 0.2)+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    ggtitle("Swa")#+

# Kil ------
ggplot(data = AO[AO$SppAlias == "Kil",], aes(x = log_ngkg_re, y = Obs_logngkg,
                      color = temperature))+
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
    facet_grid(.~system)+
    geom_encircle(aes(fill = SppAlias, color = SppAlias),
                  s_shape = 1, expand = 0, color = "black", fill = "grey",
                  alpha = 0.3,  show.legend = F)+
    geom_point(pch = 19, size = 0.2)+
    scale_color_viridis_c()+
    theme_classic()+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    ggtitle("Kil")#+

# Dar ------
ggplot(data = AO[AO$SppAlias == "Dar",], aes(x = log_ngkg_re, y = Obs_logngkg,
                      color = temperature))+
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
    facet_grid(.~system)+
    scale_color_viridis_c()+
    theme_classic()+
    geom_encircle(aes(fill = SppAlias, color = SppAlias),
                  s_shape = 1, expand = 0, color = "black", fill = "grey",
                  alpha = 0.3,  show.legend = F)+
    geom_point(pch = 19, size = 0.2)+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    ggtitle("Dar")#+

# Dac ------
ggplot(data = AO[AO$SppAlias == "Dac",], aes(x = log_ngkg_re, y = Obs_logngkg,
                      color = temperature))+
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
    facet_grid(.~system)+
    geom_encircle(aes(fill = SppAlias, color = SppAlias),
                  s_shape = 1, expand = 0, color = "black", fill = "grey",
                  alpha = 0.3,  show.legend = F)+
    geom_point(pch = 19, size = 0.2)+
    scale_color_viridis_c()+
    theme_classic()+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    ggtitle("Dac")#+


# Chu ------
ggplot(data = AO[AO$SppAlias == "Chu",], aes(x = log_ngkg_re, y = Obs_logngkg,
                      color = temperature))+
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
    facet_grid(.~system)+
    scale_color_viridis_c()+
    theme_classic()+
    geom_encircle(aes(fill = SppAlias, color = SppAlias),
                  s_shape = 1, expand = 0, color = "black", fill = "grey",
                  alpha = 0.3,  show.legend = F)+
    geom_point(pch = 19, size = 0.2)+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    ggtitle("Chu")#+

# Mad ------
ggplot(data = AO[AO$SppAlias == "Mad",], aes(x = log_ngkg_re, y = Obs_logngkg,
                      color = temperature))+
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
    facet_grid(.~system)+
    geom_encircle(aes(fill = SppAlias, color = SppAlias),
                  s_shape = 1, expand = 0, color = "black", fill = "grey",
                  alpha = 0.3,  show.legend = F)+
    geom_point(pch = 19, size = 0.2)+
    scale_color_viridis_c()+
    theme_classic()+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    ggtitle("Mad")#+

# Min ------
ggplot(data = AO[AO$SppAlias == "Min",], aes(x = log_ngkg_re, y = Obs_logngkg,
                      color = temperature))+
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
    facet_grid(.~system)+
    scale_color_viridis_c()+
    theme_classic()+
    geom_encircle(aes(fill = SppAlias, color = SppAlias),
                  s_shape = 1, expand = 0, color = "black", fill = "grey",
                  alpha = 0.3,  show.legend = F)+
    geom_point(pch = 19, size = 0.2)+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    ggtitle("Min")#+

# Pry ------
ggplot(data = AO[AO$SppAlias == "Pry",], aes(x = log_ngkg_re, y = Obs_logngkg,
                      color = temperature))+
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
    facet_grid(.~system)+
    geom_encircle(aes(fill = SppAlias, color = SppAlias),
                  s_shape = 1, expand = 0, color = "black", fill = "grey",
                  alpha = 0.3,  show.legend = F)+
    geom_point(pch = 19, size = 0.2)+
    scale_color_viridis_c()+
    theme_classic()+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')+
    ggtitle("Pry")#+

# the emergence of mean and error rates with iterations. 

for(i in 1:1000){
  
  A.0<-AO[AO$iteration <= i,] %>% 
  dplyr::group_by(system, SppAlias, PFAA) %>% 
  summarise(mean = mean(log_ngkg_re), 
            sd = sd(log_ngkg_re), 
            min = min(log_ngkg_re),
            max = max(log_ngkg_re)) %>% 
  mutate(n_iter = i) %>% 
  as.data.frame()
  
  if(i == 1){
    A.0.grow<-A.0
  }else{
    A.0.grow<-rbind(A.0.grow, A.0)
  }
}


A.3 <- AO %>%
        group_by(SppAlias, PFAA, system) %>%
        summarise(mean_val = median(Obs_logngkg, na.rm = TRUE),
                  mean_errUP = max(Obs_logngkg, na.rm = TRUE),
                  mean_errLO = min(Obs_logngkg, na.rm = TRUE))
  
  
ggplot(data = A.0.grow[A.0.grow$system == "JBA Summer" & A.0.grow$n_iter<100 &
                         !A.0.grow$SppAlias == "Phy",],
       aes(x = n_iter, y = mean,
           color = PFAA,
           group = interaction(PFAA, SppAlias)))+
  facet_grid(PFAA~SppAlias, scales="free")+
  geom_ribbon(aes(ymin = min,
                  ymax = max, color = NULL, fill = PFAA), alpha = 0.4) +
  geom_ribbon(aes(ymin = mean - sd,
                  ymax = mean + sd, color = NULL, fill = PFAA), alpha = 1) + 
  geom_line(color = "black")+
  geom_errorbar(data = A.3[A.3$system == "JBA Summer" & !A.3$SppAlias == "Phy",],
          mapping = aes(y = mean_val, x = -5,
                        ymax = mean_errUP,
                        ymin = mean_errLO,
                        color = PFAA, 
           group = interaction(PFAA, SppAlias)))+
  geom_point(data = A.3[A.3$system == "JBA Summer" & !A.3$SppAlias == "Phy",],
          mapping = aes(y = mean_val, x = -5,
                        ymax = mean_errUP,
                        ymin = mean_errLO,
                        fill = PFAA, 
           group = interaction(PFAA, SppAlias)), pch=21, color = "black")+
  theme_classic()
# add sample sizes 



ggplot(data = A.0.grow[A.0.grow$system == "JBA Spring" & A.0.grow$n_iter<100,], aes(x = n_iter, y = mean,
                            color = PFAA,
                            group = interaction(PFAA, SppAlias)))+
  facet_grid(PFAA~SppAlias, scales="free")+
  geom_ribbon(aes(ymin = min,
                  ymax = max, color = NULL, fill = PFAA), alpha = 0.4) +
  geom_ribbon(aes(ymin = mean - sd,
                  ymax = mean + sd, color = NULL, fill = PFAA), alpha = 1) + 
  geom_line(color = "black")+
  theme_classic()

ggplot(data = A.0.grow[A.0.grow$system == "WG Fall" & A.0.grow$n_iter<100,], aes(x = n_iter, y = mean,
                            color = PFAA,
                            group = interaction(PFAA, SppAlias)))+
  facet_grid(PFAA~SppAlias, scales="free")+
  geom_ribbon(aes(ymin = min,
                  ymax = max, color = NULL, fill = PFAA), alpha = 0.4) +
  geom_ribbon(aes(ymin = mean - sd,
                  ymax = mean + sd, color = NULL, fill = PFAA), alpha = 1) + 
  geom_line(color = "black")+
  theme_classic()

# sediment water ratio
AO$wto_s_pfas_ratio <- AO$C_WTO/ AO$C_s

ggformat( plot = ggplot(AO[!c(AO$SppAlias == "Phy"),],
       aes(y=log_ngkg_re, x = wto_s_pfas_ratio,
           group = interaction( PFAA, spp_PFAA),
           color = PFAA, size = C_s))+
  # geom_smooth(method = "lm", se = F)+
  # geom_point( pch=21)+
  geom_line(linewidth = 0.1)+
  facet_grid(system~SppAlias, scales = "free")+
  scale_color_manual(values = PFAA_cols)+
  scale_linetype_manual(values = c("dashed", "dotted", "solid"))+
  theme_bw(),
 x_title = ("PFAA water:sediment ratio"),
         y_title = ("Tissue PFAA (ng/kg)"))

ggformat( plot = ggplot(AO[!c(AO$SppAlias == "Phy") & AO$PFAA == "PFOS",],
       aes(y=Dietary_uptake, x = wto_s_pfas_ratio,
           group = interaction( PFAA, spp_PFAA),
           color = PFAA))+
  # geom_smooth(method = "lm", se = F)+
  # geom_point( pch=".")+
  geom_line(linewidth = 0.1)+
  facet_grid(system~PFAA, scales = "free")+
  scale_color_manual(values = PFAA_cols)+
  scale_linetype_manual(values = c("dashed", "dotted", "solid"))+
  theme_bw(),
 x_title = ("PFAA water:sediment ratio"),
         y_title = ("Dietary uptake"))

ggformat( plot = ggplot(AO[!c(AO$SppAlias == "Phy") & AO$PFAA == "PFOA",],
       aes(y=Dietary_uptake, x = wto_s_pfas_ratio,
           group = interaction( PFAA, spp_PFAA),
           color = PFAA))+
  # geom_smooth(method = "lm", se = F)+
  # geom_point( pch=".")+
  geom_line(linewidth = 0.1)+
  facet_grid(system~PFAA, scales = "free")+
  scale_color_manual(values = PFAA_cols)+
  scale_linetype_manual(values = c("dashed", "dotted", "solid"))+
  theme_bw(),
 x_title = ("PFAA water:sediment ratio"),
         y_title = ("Dietary uptake"))

ggformat( plot = ggplot(AO[!c(AO$SppAlias == "Phy") & AO$PFAA == "PFHxS",],
       aes(y=Dietary_uptake, x = wto_s_pfas_ratio,
           group = interaction( PFAA, spp_PFAA),
           color = PFAA))+
  # geom_smooth(method = "lm", se = F)+
  # geom_point( pch=".")+
  geom_line(linewidth = 0.1)+
  facet_grid(system~PFAA, scales = "free")+
  scale_color_manual(values = PFAA_cols)+
  scale_linetype_manual(values = c("dashed", "dotted", "solid"))+
  theme_bw(),
 x_title = ("PFAA water:sediment ratio"),
         y_title = ("Dietary uptake"))

# base plot with mean plots and ranges as error bars
A2<-AO %>% 
  dplyr::group_by(system, SppAlias, PFAA) %>% 
  summarise(mean_pred = mean(log_ngkg_re), 
            sd_pred = sd(log_ngkg_re), 
            min_pred = min(log_ngkg_re),
            max_pred = max(log_ngkg_re), 
            mean_obs = mean(Obs_logngkg), 
            sd_obs = sd(Obs_logngkg), 
            min_obs = min(Obs_logngkg),
            max_obs = max(Obs_logngkg)) %>% 
  as.data.frame()

ggplot(data = A2, aes(x = mean_pred, y = mean_obs, color = PFAA))+
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
    facet_grid(.~system)+
    # geom_errorbar(aes(ymin = min_obs, ymax = max_obs), linewidth = 0.5) +
    # geom_errorbarh(aes(xmin = min_pred, xmax = max_pred), linewidth = 0.5) + 
    # geom_point(pch = 19, size = 2)+
    theme_classic()+
    geom_text(aes(label = SppAlias), size = 2)+
    scale_color_manual(values = PFAA_cols)+
    scale_fill_manual(values = PFAA_cols)+
    xlab('Modeled log PFAA (ng/kg)')+
    ylab('Observed log PFAA (ng/kg)')#+
    # coord_cartesian(ylim = c(2,6.5), xlim = c(2,6.5))+
    # theme(legend.position = "none")
# p1_fish




