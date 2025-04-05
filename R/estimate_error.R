
#' Title: Model bias and simple model error estimates 
#' Use: 
#' 
#' @param MetricData 
#' @param DataID 
#' @param type two options: "ArnotGobas2004" (default) and "McLeodetal2016"
#' @param verbose logical to print messages or not. 
#'
#' @return
#' @export
#'
#' 
#' 
estimate_error <- function(MetricData,
                           DataID,
                           type = "ArnotGobas2004",
                           verbose = TRUE){
    
  MetricData <- MetricData[c(!is.na(MetricData$Obs_ngg) | !is.nan(MetricData$Obs_ngg)), ] ## drop values with ND observations not shown on plots

  # Model Bias (Arnot and Gobas 2004)  -------------
  if(type == "ArnotGobas2004"){
    ## Calculate Model Bias by PFAA
    n.PFAA <- length(unique(MetricData$PFAA))
    m.SppAlias<- length(unique(MetricData$SppAlias))
    
    ## Calculate MB for individual species 
    # ********************************************************
    MetricData$MB_Tissue <- log10((10^MetricData[, 'logC_B_ngkg']) / (10^MetricData[,'Obs_logngkg']))
    MetricData$MB_Tissue.n <- MetricData$MB_Tissue/n.PFAA
    
    MBj<-MetricData %>% 
      dplyr::group_by(SppAlias) %>%
      summarize(MBj = 10^(sum(MB_Tissue.n))) 
    
    ## Calculate MB for individual PFAA
    # ********************************************************
    MetricData$MB_Tissue <- log10((10^MetricData[, 'logC_B_ngkg']) / (10^MetricData[,'Obs_logngkg']))
    MetricData$MB_Tissue.n2 <- MetricData$MB_Tissue/m.SppAlias 
    
    MBi<-MetricData %>% 
      dplyr::group_by(PFAA) %>%
      summarize(MBi = 10^(sum(MB_Tissue.n2, na.rm = TRUE))) 
    
    ## The entire Model Bias
    # ******************************************************** 
    MB<-10^(sum(MBi$MBi/m.SppAlias)) 
    
  } else if (type == "McLeodetal2016"){
    # Model Bias (McLeod et al 2016; DOI: 10.1021/acs.est.6b03169) -------------
    # model MB_avg bias used to compare the overall predictions 
    # model MB_cv (coefficient of variation) used to compare effectiveness in each model 
    # for predicting variation
  
    MB_avg <- (log(mean(MetricData[,'Obs_ngg'], na.rm = TRUE))/# multiply to make it ngkg
                 log(mean(MetricData[,'C_B_ngg'], na.rm = TRUE)))
    
    cv_obs <- sd(MetricData[,'Obs_ngg'], na.rm = TRUE)/
      mean(MetricData[,'Obs_ngg'], na.rm = TRUE) # multiply to make it ngkg
    
    cv_model <- sd(MetricData[,'C_B_ngg'], na.rm = TRUE)/
      mean(MetricData[,'C_B_ngg'], na.rm = TRUE)
    
    MB_CV <-  (log(cv_obs)/log(cv_model)) 
      
    # Model Bias (McLeod et al 2016 - per species)
    MB_avg_species<- MetricData %>% 
      dplyr::group_by(SppAlias) %>%
      summarise(log(mean(Obs_ngg, na.rm = TRUE))/log(mean(C_B_ngg, na.rm = TRUE)))
    
    # Model Bias (McLeod et al 2016 - per species)
    MB_avg_PFAA<- MetricData %>% 
      dplyr::group_by(PFAA) %>%
      summarise(log(mean(Obs_ngg, na.rm = TRUE))/log(mean(C_B_ngg, na.rm = TRUE)))
  
  }else {
    message("No model bias estimates, 'type' must be either 'McLeodetal2016' or 'ArnotGobas2004'")
  }
  
  # % bias calculations (Kraskura and Salice) -------------

  # overall model
  MB_avg_perc <- (mean(MetricData[,'Obs_ngg'], na.rm = TRUE)-
                    (mean(MetricData[,'C_B_ngg'], na.rm = TRUE)))/
                 mean(MetricData[,'Obs_ngg'], na.rm = TRUE)  * 100

  # Model Bias (McLeod et al 2016 - per species)
  MB_avg_perc_species<- MetricData %>% 
    dplyr::group_by(SppAlias) %>%
    summarise((mean(Obs_ngg, na.rm = TRUE)/mean(C_B_ngg, na.rm = TRUE))/
                mean(Obs_ngg, na.rm = TRUE) * 100)
  
  # Model Bias (McLeod et al 2016 - per species)
  MB_avg_perc_PFAA<- MetricData %>% 
    dplyr::group_by(PFAA) %>%
    summarise((mean(Obs_ngg, na.rm = TRUE)/mean(C_B_ngg, na.rm = TRUE))/
                mean(Obs_ngg, na.rm = TRUE) * 100)
  
  # RSS and TSS, R2 -----
  MetricData <- MetricData[!c(is.infinite(MetricData$Obs_logngkg) | is.na(MetricData$Obs_logngkg)),]
  
  rss <- sum((MetricData$logC_B_ngkg - MetricData$Obs_logngkg) ^ 2)  ## residual sum of squares
  tss <- sum((MetricData$Obs_logngkg - mean(MetricData$Obs_logngkg)) ^ 2, na.rm = TRUE)  ## total sum of squares
  rsq <- 1 - rss/tss
  
  # r2 by PFAA (not used)
  r2_pfaa <- MetricData %>%
    dplyr::group_by(PFAA) %>%
    summarise(mean_Obs_logngkg = mean(Obs_logngkg),
              RSS = sum((logC_B_ngkg - Obs_logngkg) ^ 2), 
              TSS = sum((logC_B_ngkg - mean_Obs_logngkg) ^ 2),
              r2 = 1 - (RSS/TSS)) %>% 
    as.data.frame()
  
  # r2 species (not used)
  r2_species <- MetricData %>%
    dplyr::group_by(SppAlias) %>%
    summarise(mean_Obs_logngkg = mean(Obs_logngkg),
              RSS = sum((logC_B_ngkg - Obs_logngkg) ^ 2), 
              TSS = sum((logC_B_ngkg - mean_Obs_logngkg) ^ 2),
              r2 = 1 - (RSS/TSS)) %>% 
    as.data.frame()

  colnames(MB_avg_species) <- c("SppAlias", "MBavg")
  colnames(MB_avg_PFAA) <- c("PFAA", "MBavg")


    message(paste("Using", type, "approach of estimating model bias"))

    if(type == "ArnotGobas2004"){
      
      if(verbose){message("returns: MB,  r2, MBi, MBj,
          MB_avg_perc, MB_avg_perc_species, MB_avg_perc_PFAA")}
      
      return(list(MB_avg, rsq, MBi, MBj,
              MB_avg_perc, MB_avg_perc_species, MB_avg_perc_PFAA))
      
    }else if(type == "McLeodetal2016"){
      
      if(verbose){message("returns: MBavg, MBCV, r2, MB_avg_species, MB_avg_PFAA,
          MB_avg_perc, MB_avg_perc_species, MB_avg_perc_PFAA")}
      
      return(list(MB_avg, MB_CV, rsq, MB_avg_species, MB_avg_PFAA,
              MB_avg_perc, MB_avg_perc_species, MB_avg_perc_PFAA))
      
    }
}
