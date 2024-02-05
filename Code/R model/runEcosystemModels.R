# loop without RE and with RE: output MEDIAN (0) RESULTS  
PFAA_loop <- function(settings,
                     inputFiles_list,
                     parameterList,
                     PFAA_List){

    for (i in 1:length(PFAA_List)){
        PFAA<-PFAA_List[[i]]
        # TissueConc <- BioaccumulationModel(PFAA = PFAA, settings_dict = settings, **inputFiles_list)
        TissueConc_re <- BioaccumulationModel(
                          PFAA = PFAA,
                          numSpecies = inputFiles_list$numSpecies,
                          oceanData = inputFiles_list$oceanData,
                          chemicalData = inputFiles_list$chemicalData,
                          chemicalParams = inputFiles_list$chemicalParams,
                          organismData = inputFiles_list$organismData,
                          foodWebData = inputFiles_list$foodWebData,
                          settings = settings,
                          dietData = inputFiles_list$dietData,
                          kRTable = inputFiles_list$kRTable)
          temp_table_re <- TissueConc_re[parameterList]
          temp_table_re$PFAA <- PFAA
          if(i == 1){
            ResultTable_re <- temp_table_re
          } else {
            ResultTable_re<-rbind(ResultTable_re, temp_table_re)
          }
    }
  
    ResultTable_re$SppAlias <- rownames(ResultTable_re)
    rownames(ResultTable_re) <- 1:nrow(ResultTable_re)
    names(ResultTable_re)[names(ResultTable_re) == 'C_B'] <- 'C_B_re' # ng/kg
    ResultTable_re$log_ngkg_re <- log10(ResultTable_re$C_B_re)
    ResultTable_re$log_BAF_re <- log10(ResultTable_re$C_B_re / (ResultTable_re$C_WTO*1000))
    
    ResultTable_re$SppAlias <- substr(ResultTable_re$SppAlias, start = 1, stop =3)
    ResultTable_re<-ResultTable_re[, c("SppAlias", "PFAA" ,"C_B_re", parameterList[2:length(parameterList)],  "log_BAF_re", "log_ngkg_re")]
  
    AllResults<-ResultTable_re
    return(AllResults)
  
}

runEcosystemModel <- function(settings,
                              inputFiles_list, 
                              PFAA_List,
                              kRTable = NULL, 
                              dietData = NULL, 
                              parameterList = NULL, 
                              min_diet_WTO_C_s = FALSE, 
                              max_diet_WTO_C_s = FALSE, 
                              min_dietData = NULL, 
                              max_dietData = NULL, 
                              observedData_BCF_BMF = NULL,
                              env_temperature = NULL, # only works for food web full model
                              env_pH = NULL,
                              env_DO = NULL,
                              env_Cs = NULL, # as many PFAS as in the dataset
                              env_WTO = NULL, # as many PFAS as in the dataset
                              org_mass = NULL, # as many fish as in the dataset 
                              RunID = "default",
                              random_temp = FALSE,
                              random_DO = FALSE
                              ){ # for BCF_calc or BMF_calc = TRUE){ 
  
  # food web calc = TRUE ----


  inputFiles_list$dietData <- dietData
  inputFiles_list$kRTable <- kRTable
  
  if(!is.null(org_mass)){
    for(i in 1:length(org_mass)){
      inputFiles_list$organismData[1,i+1] <- org_mass[i]
    }
  }

  if(!is.null(env_temperature)){
    inputFiles_list$oceanData["T", 1] <- env_temperature
  }
  if(!is.null(env_pH)){
    inputFiles_list$oceanData["pH", 1] <- env_pH
  }
  if(!is.null(env_DO)){
    inputFiles_list$oceanData["C_OX", 1] <- env_DO
  }
  
  if(!is.null(env_WTO)){
    inputFiles_list$chemicalData["C_WTO", 2:(length(env_WTO)+1)] <- env_WTO
  }
  if(!is.null(env_Cs)){
    inputFiles_list$chemicalData["C_s", 2:(length(env_Cs)+1)] <- env_Cs
  }

  run0<-PFAA_loop(settings,
                 inputFiles_list,
                 parameterList,
                 PFAA_List)
  run0.full<-as.data.frame(run0)

  if(!is.null(max_dietData)){
    
    if(max_diet_WTO_C_s){ # max diet and max observed water and sediment samples
      chemicalData <- inputFiles_list$chemicalData
      rownames(chemicalData) <- c("C_WTO_median", "C_WTO", "C_WTO_min",
                                  "C_s_median", "C_s", "C_s_min", "Phi_exp")
      inputFiles_list$chemicalData <- chemicalData
    }
    
    inputFiles_list$dietData <- max_dietData
    run.MAX<-PFAA_loop(settings,
                       inputFiles_list,
                       parameterList,
                       PFAA_List)
    run0.MAXfull<-as.data.frame(run.MAX)
    
    names(run0.MAXfull)[names(run0.MAXfull) == 'log_ngkg_re'] <- 'log_ngkg_re_max'
    names(run0.MAXfull)[names(run0.MAXfull) == 'C_s'] <- 'C_s_max'
    names(run0.MAXfull)[names(run0.MAXfull) == 'C_WTO'] <- 'C_WTO_max'
    names(run0.MAXfull)[names(run0.MAXfull) == 'C_B_re'] <- 'C_B_re_max'
    names(run0.MAXfull)[names(run0.MAXfull) == 'Gill_uptake'] <- 'Gill_uptake_max'
    names(run0.MAXfull)[names(run0.MAXfull) == 'Dietary_uptake'] <- 'Dietary_uptake_max'
    names(run0.MAXfull)[names(run0.MAXfull) == "log_BAF_re"] <- 'log_BAF_re_max'
  }
  
  if(!is.null(min_dietData)){
    
    if(min_diet_WTO_C_s){
      chemicalData <- inputFiles_list$chemicalData
      rownames(chemicalData) <- c("C_WTO_median", "C_WTO_max", "C_WTO", 
                                  "C_s_median", "C_s_max", "C_s", "Phi_exp")
      inputFiles_list$chemicalData <- chemicalData        
    }

    inputFiles_list$dietData <- min_dietData

    run.MIN<-PFAA_loop(settings,
                       inputFiles_list,
                       parameterList,
                       PFAA_List)
    run0.MINfull<-as.data.frame(run.MIN)
  
    names(run0.MINfull)[names(run0.MINfull) == 'log_ngkg_re'] <- 'log_ngkg_re_min'
    names(run0.MINfull)[names(run0.MINfull) == 'C_s'] <- 'C_s_min'
    names(run0.MINfull)[names(run0.MINfull) == 'C_WTO'] <- 'C_WTO_min'
    names(run0.MINfull)[names(run0.MINfull) == 'C_B_re'] <- 'C_B_re_min'
    names(run0.MINfull)[names(run0.MINfull) == 'Gill_uptake'] <- 'Gill_uptake_min'
    names(run0.MINfull)[names(run0.MINfull) == 'Dietary_uptake'] <- 'Dietary_uptake_min'
    names(run0.MINfull)[names(run0.MINfull) == "log_BAF_re"] <- 'log_BAF_re_min'
  }
  
  # reset the median, mean and max values and names in the dataset
  if(!is.null(min_dietData) && !is.null(max_dietData)){
    if(max_diet_WTO_C_s & min_diet_WTO_C_s){
      # reset chemicalData
      chemicalData <- inputFiles_list$chemicalData
      rownames(chemicalData) <- c("C_WTO", "C_WTO_max", "C_WTO_min",
                                  "C_s", "C_s_max", "C_s_min", "Phi_exp")
      inputFiles_list$chemicalData <- chemicalData        
    }
    
    run0.full<-run0.full[ , !(names(run0.full) %in% c("Total Elimination.1"))]
    run0.MAXfull<-run0.MAXfull[ , !(names(run0.MAXfull) %in% c("Total Elimination.1"))]
    run0.MINfull<-run0.MINfull[ , !(names(run0.MINfull) %in% c("Total Elimination.1"))]

    
    run0.MAX.MIN <- merge(run0.MAXfull, run0.MINfull, by=c('SppAlias','PFAA',
                                                           "k1", "k2", "ke", 
                                                           "kg", "kd", "RMR",
                                                           "G_V", "G_D",
                                                           "Total Elimination"))
    Modeled <- merge(run0.full, run0.MAX.MIN, by=c('SppAlias','PFAA',
                                                   "k1", "k2", "ke", 
                                                   "kg", "kd", "RMR",
                                                   "G_V", "G_D",
                                                   "Total Elimination"))

    Modeled$logngkg_re_Up <- Modeled$log_ngkg_re_max - Modeled$log_ngkg_re
    Modeled$logngkg_re_Lo <- Modeled$log_ngkg_re - Modeled$log_ngkg_re_min
      
    # Format observational (median dietData) data 
    # dietData == the same as the median data
    ObsAvg <- dietData[!(row.names(dietData) %in% c('Phy')),]
    ObsAvg$SppAlias <- rownames(ObsAvg)
    ObsAvg <- gather(ObsAvg, PFAA, Obs_ngg, PFHxS:PFUA)
    ObsAvg$Obs_logngkg <- log10(ObsAvg$Obs_ngg * 1000)
    
    ObsMax <- max_dietData[!(row.names(max_dietData) %in% c('Phy')),]
    ObsMax$SppAlias <- rownames(ObsMax)
    ObsMax <- gather(ObsMax, PFAA, Obs_ngg_max, PFHxS:PFUA)
    ObsMax$Obs_logngkg_max <- log10(ObsMax$Obs_ngg_max * 1000)

    ObsMin <- min_dietData[!(row.names(min_dietData) %in% c('Phy')),]
    ObsMin$SppAlias <- rownames(ObsMin)
    ObsMin <- gather(ObsMin, PFAA, Obs_ngg_min, PFHxS:PFUA)
    ObsMin$Obs_logngkg_min <- log10(ObsMin$Obs_ngg_min * 1000)
    
    Observed <- merge(merge(ObsAvg, ObsMax, by=c('SppAlias','PFAA'), all = TRUE), ObsMin, by=c('SppAlias','PFAA'), all = TRUE)
    Observed$Obs_logngkg_Up <- Observed$Obs_logngkg_max - Observed$Obs_logngkg # median to the max
    Observed$Obs_logngkg_Lo <- Observed$Obs_logngkg - Observed$Obs_logngkg_min # min to the median
    AllData <- merge(Modeled, Observed, all = TRUE, by=c("SppAlias",'PFAA'))
  
  }else{ 
  # if diet data min and max do not exist
    run0.full$log_ngkg_re_max<- NA
    run0.full$log_ngkg_re_min<- NA
    run0.full$logngkg_re_Up<- NA
    run0.full$logngkg_re_Lo<- NA
    Modeled<-run0.full
    
    # Format observational (median dietData) data 
    # dietData == the same as the median data
    ObsAvg <- dietData[!(row.names(dietData) %in% c('Phy')),]
    ObsAvg$SppAlias <- rownames(ObsAvg)
    ObsAvg <- reshape::melt(ObsAvg, id_vars=c('SppAlias'), variable.name = 'PFAA', value.name ='Obs_ngg' )
    colnames(ObsAvg) <- c("SppAlias", "PFAA", "Obs_ngg")
    ObsAvg$Obs_logngkg <- log10(ObsAvg$Obs_ngg * 1000) # ng/kg conversion
    AllData <- merge(Modeled, ObsAvg, all = TRUE, by = c("SppAlias",'PFAA'))
    
    AllData$Obs_ngg_max <- NA
    AllData$Obs_logngkg_max <- NA
    AllData$Obs_ngg_min <- NA
    AllData$Obs_logngkg_min <- NA
    AllData$Obs_logngkg_Up <- NA
    AllData$Obs_logngkg_Lo <- NA
  }
  
  AllData$temperature<- inputFiles_list$oceanData["T", 1] # only works for food web full model
  AllData$pH<- inputFiles_list$oceanData["pH", 1] 
  AllData$env_DO <- inputFiles_list$oceanData["C_OX", 1]
  AllData$runID <- RunID
  species_data<-as.data.frame(t(inputFiles_list$organismData[,2:inputFiles_list$numSpecies])[,1])
  species_data$SppAlias<-rownames(species_data)
  colnames(species_data)<-c("WB", "SppAlias")
  AllData<-merge(AllData, species_data, by  = "SppAlias", all.x = T)
  
  
  # returns ----
  return(AllData)

  
}
