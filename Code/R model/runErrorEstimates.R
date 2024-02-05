
# whole model MB. 

estimate_error <- function(MetricData, DataID){
  ## Calculate Model Bias by PFAA
  MetricData <- MetricData[c(!is.na(MetricData$Obs_ngg) | !is.nan(MetricData$Obs_ngg)), ] ## drop values with ND observations not shown on plots
  n.PFAA <- length(unique(MetricData$PFAA))
  m.SppAlias<- length(unique(MetricData$SppAlias))
  
  
  ## Calculate MB for individual species 
  # ********************************************************
  MetricData$MB_Tissue <- log10((10^MetricData[, 'log_ngkg_re']) / (10^MetricData[,'Obs_logngkg']))
  MetricData$MB_Tissue.n <- MetricData$MB_Tissue/n.PFAA
  
  MBj<-MetricData %>% 
    dplyr::group_by(SppAlias) %>%
    summarize(MBj = 10^(sum(MB_Tissue.n))) 
  
  ## Calculate MB for individual PFAA
  # ********************************************************
  MetricData$MB_Tissue <- log10((10^MetricData[, 'log_ngkg_re']) / (10^MetricData[,'Obs_logngkg']))
  MetricData$MB_Tissue.n2 <- MetricData$MB_Tissue/m.SppAlias
  
  MBi<-MetricData %>% 
    dplyr::group_by(PFAA) %>%
    summarize(MBi = 10^(sum(MB_Tissue.n2, na.rm = TRUE))) 
  
  ## The entire Model Bias
  # ******************************************************** 
  MB<-10^(sum(MBi$MBi/m.SppAlias))
  
  # RSS and TSS, R2 -----
  MetricData <- MetricData[!c(is.infinite(MetricData$Obs_logngkg) | is.na(MetricData$Obs_logngkg)),]
  
  rss <- sum((MetricData$log_ngkg_re - MetricData$Obs_logngkg) ^ 2)  ## residual sum of squares
  tss <- sum((MetricData$Obs_logngkg - mean(MetricData$Obs_logngkg)) ^ 2, na.rm = TRUE)  ## total sum of squares
  rsq <- 1 - rss/tss
  
  # r2 by PFAA
  r2_pfaa <- MetricData %>%
    dplyr::group_by(PFAA) %>%
    summarise(RSS = sum((log_ngkg_re - Obs_logngkg) ^ 2)) %>% 
    as.data.frame()
  
  # r2 species
  r2_species <- MetricData %>%
    dplyr::group_by(SppAlias) %>%
    summarise(RSS = sum((log_ngkg_re - Obs_logngkg) ^ 2)) %>% 
    as.data.frame()

  # colnames(ModelMetricsSpp) <- c("SppAlias", paste(DataID, "MBxTissue"))

  return(list(MBi, MBj, MB, rsq,  r2_pfaa, r2_species))
}
