#' Title: 
#' Description: Runs Bioaccumulation Model function which estimates PFAS (one compound at a time) in each species across the food web.
#' 
#' @param settings 
#' @param inputFiles_list listed input data frames (provided by the author, not calculated)
#' @param parameterList 
#' @param PFAA_List 
#' @param median_values
#' @param min_values
#' @param max_values
#'
#' @returns
#' @export
#'
#' @examples
#' 
#' 
food_web_bioaccumulation <- function(settings,
                      inputFiles_list,
                      parameterList,
                      PFAA_List,
                      median_values = FALSE,
                      min_values = FALSE,
                      max_values = FALSE){
  # mdian values
  if(median_values){
    for (i in 1:length(PFAA_List)){
        PFAA<-PFAA_List[[i]]

        TissueConc_re <- bioaccumulation_model(
                          PFAA = PFAA,
                          numSpecies = inputFiles_list$numSpecies,
                          oceanData = inputFiles_list$oceanData,
                          chemicalData = inputFiles_list$chemicalData,
                          chemicalParams = inputFiles_list$chemicalParams,
                          organismData = inputFiles_list$organismData,
                          foodWebData = inputFiles_list$foodWebData,
                          settings = settings,
                          dietData = inputFiles_list$dietData,
                          kRTable = inputFiles_list$kRTable,
                          median_values = TRUE)
          # temp_table_re <- TissueConc_re[parameterList]
          temp_table_re <- TissueConc_re
          temp_table_re$PFAA <- PFAA
          if(i == 1){
            ResultTable_re <- temp_table_re
          } else {
            ResultTable_re<-rbind(ResultTable_re, temp_table_re)
          }
    }

    # Add temperature, pH, oxygen conc (C_OX), and species masses
    ResultTable_re$temperature<- inputFiles_list$oceanData["T", 1] 
    ResultTable_re$pH<- inputFiles_list$oceanData["pH", 1]
    ResultTable_re$env_DO <- inputFiles_list$oceanData["C_OX", 1]
    
    species_data<-as.data.frame(t(inputFiles_list$organismData[,2:inputFiles_list$numSpecies])[,1])
    species_data$SppAlias<-rownames(species_data)
    colnames(species_data)<-c("WB", "SppAlias")
    
    ResultTable_re<-merge(ResultTable_re, species_data, by  = "SppAlias", all.x = T)

    # Format observational (median dietData) data
    dietData <- as.data.frame(inputFiles_list$dietData)
    ObsData <- dietData[!(row.names(dietData) %in% c('Phy', 'Plt')),]
    ObsData$SppAlias <- rownames(ObsData)
    ObsData<-(tidyr::pivot_longer(ObsData,
                        cols = c(1:length(PFAA_List)),      
                        id_vars = 'SppAlias',
                        names_to = 'PFAA',
                        values_to ='Obs_ngg'))
    # print(ObsData)
    colnames(ObsData) <- c("SppAlias", "PFAA", "Obs_ngg")
    ObsData$Obs_ngkg <- ObsData$Obs_ngg * 1000 # ng/kg conversion
    ObsData$Obs_logngkg <- log10(ObsData$Obs_ngg * 1000) # ng/kg conversion
    
    ResultTable_re <- merge(ResultTable_re, ObsData, all = TRUE, by = c("SppAlias",'PFAA'))

    # Select parameters as identified in 'parameterList'
    ResultTable_re_median<-ResultTable_re[, c("SppAlias", "PFAA" ,"C_B_ngkg","C_B_ngg","logC_B_ngkg",
                                       "Obs_ngg","Obs_logngkg", "Obs_ngkg",
                                       parameterList[1:length(parameterList)])]
    ResultTable_re_median$Spp_PFAA<-paste(ResultTable_re_median$SppAlias, ResultTable_re_median$PFAA, sep = "-")

    message("model run with median values of PFAS and environmental observations")
  
  } # median value block end
  
  # min values 
  if(min_values){
    for (i in 1:length(PFAA_List)){
        PFAA<-PFAA_List[[i]]

        TissueConc_re <- bioaccumulation_model(
                          PFAA = PFAA,
                          numSpecies = inputFiles_list$numSpecies,
                          oceanData = inputFiles_list$oceanData,
                          chemicalData = inputFiles_list$chemicalData,
                          chemicalParams = inputFiles_list$chemicalParams,
                          organismData = inputFiles_list$organismData,
                          foodWebData = inputFiles_list$foodWebData,
                          settings = settings,
                          dietData = inputFiles_list$min_dietData,
                          kRTable = inputFiles_list$kRTable,
                          min_values = TRUE)
          # temp_table_re <- TissueConc_re[parameterList]
          temp_table_re <- TissueConc_re
          temp_table_re$PFAA <- PFAA
          if(i == 1){
            ResultTable_re <- temp_table_re
          } else {
            ResultTable_re<-rbind(ResultTable_re, temp_table_re)
          }
    }

    # Add temperature, pH, oxygen conc (C_OX), and species masses
    ResultTable_re$temperature<- inputFiles_list$oceanData["T", 1] 
    ResultTable_re$pH<- inputFiles_list$oceanData["pH", 1]
    ResultTable_re$env_DO <- inputFiles_list$oceanData["C_OX", 1]
    
    species_data<-as.data.frame(t(inputFiles_list$organismData[,2:inputFiles_list$numSpecies])[,1])
    species_data$SppAlias<-rownames(species_data)
    colnames(species_data)<-c("WB", "SppAlias")
    
    ResultTable_re<-merge(ResultTable_re, species_data, by  = "SppAlias", all.x = T)

    # Format observational (median dietData) data
    dietData <- as.data.frame(inputFiles_list$min_dietData)
    ObsData <- dietData[!(row.names(dietData) %in% c('Phy', 'Plt')),]
    ObsData$SppAlias <- rownames(ObsData)
    ObsData<-(tidyr::pivot_longer(ObsData,
                        cols = c(1:length(PFAA_List)),      
                        id_vars = 'SppAlias',
                        names_to = 'PFAA',
                        values_to ='Obs_ngg'))
    # print(ObsData)
    colnames(ObsData) <- c("SppAlias", "PFAA", "Obs_ngg")
    ObsData$Obs_ngkg <- ObsData$Obs_ngg * 1000 # ng/kg conversion
    ObsData$Obs_logngkg <- log10(ObsData$Obs_ngg * 1000) # ng/kg conversion
    
    ResultTable_re <- merge(ResultTable_re, ObsData, all = TRUE, by = c("SppAlias",'PFAA'))

    # Select parameters as identified in 'parameterList'
    ResultTable_re_min<-ResultTable_re[, c("SppAlias", "PFAA" ,"C_B_ngkg","C_B_ngg","logC_B_ngkg",
                                       "Obs_ngg","Obs_logngkg", "Obs_ngkg",'C_WTO','C_s')]
    ResultTable_re_min$Spp_PFAA<-paste(ResultTable_re_min$SppAlias, ResultTable_re_min$PFAA, sep = "-")
    ResultTable_re_min <- ResultTable_re_min %>% 
      dplyr::rename_at(vars(-c("Spp_PFAA")), function(x) paste0(x, "_min"))

    message("model run with min values of PFAS and environmental observations")
  
  } # min value block end

    # maxvalues 
  if(max_values){
    for (i in 1:length(PFAA_List)){
        PFAA<-PFAA_List[[i]]

        TissueConc_re <- bioaccumulation_model(
                          PFAA = PFAA,
                          numSpecies = inputFiles_list$numSpecies,
                          oceanData = inputFiles_list$oceanData,
                          chemicalData = inputFiles_list$chemicalData,
                          chemicalParams = inputFiles_list$chemicalParams,
                          organismData = inputFiles_list$organismData,
                          foodWebData = inputFiles_list$foodWebData,
                          settings = settings,
                          dietData = inputFiles_list$max_dietData,
                          kRTable = inputFiles_list$kRTable,
                          max_values = TRUE)
          # temp_table_re <- TissueConc_re[parameterList]
          temp_table_re <- TissueConc_re
          temp_table_re$PFAA <- PFAA
          if(i == 1){
            ResultTable_re <- temp_table_re
          } else {
            ResultTable_re<-rbind(ResultTable_re, temp_table_re)
          }
    }

    # Add temperature, pH, oxygen conc (C_OX), and species masses
    ResultTable_re$temperature<- inputFiles_list$oceanData["T", 1] 
    ResultTable_re$pH<- inputFiles_list$oceanData["pH", 1]
    ResultTable_re$env_DO <- inputFiles_list$oceanData["C_OX", 1]
    
    species_data<-as.data.frame(t(inputFiles_list$organismData[,2:inputFiles_list$numSpecies])[,1])
    species_data$SppAlias<-rownames(species_data)
    colnames(species_data)<-c("WB", "SppAlias")
    
    ResultTable_re<-merge(ResultTable_re, species_data, by  = "SppAlias", all.x = T)

    # Format observational (median dietData) data
    dietData <- as.data.frame(inputFiles_list$max_dietData)
    ObsData <- dietData[!(row.names(dietData) %in% c('Phy', 'Plt')),]
    ObsData$SppAlias <- rownames(ObsData)
    ObsData<-(tidyr::pivot_longer(ObsData,
                        cols = c(1:length(PFAA_List)),      
                        id_vars = 'SppAlias',
                        names_to = 'PFAA',
                        values_to ='Obs_ngg'))
    # print(ObsData)
    colnames(ObsData) <- c("SppAlias", "PFAA", "Obs_ngg")
    ObsData$Obs_ngkg <- ObsData$Obs_ngg * 1000 # ng/kg conversion
    ObsData$Obs_logngkg <- log10(ObsData$Obs_ngg * 1000) # ng/kg conversion
    
    ResultTable_re <- merge(ResultTable_re, ObsData, all = TRUE, by = c("SppAlias",'PFAA'))

    # Select parameters as identified in 'parameterList'
    ResultTable_re_max<-ResultTable_re[, c("SppAlias", "PFAA" ,"C_B_ngkg","C_B_ngg","logC_B_ngkg",
                                       "Obs_ngg","Obs_logngkg", "Obs_ngkg",'C_WTO','C_s')]
    ResultTable_re_max$Spp_PFAA<-paste(ResultTable_re_max$SppAlias, ResultTable_re_max$PFAA, sep = "-")
    
    ResultTable_re_max <- ResultTable_re_max %>% 
      dplyr::rename_at(vars(-c("Spp_PFAA")), function(x) paste0(x, "_max"))
    
    
    message("model run with max values of PFAS and environmental observations")
  
  } # max value block end
  
  if(median_values & !c(min_values | max_values)){
    ResultTable_re <- ResultTable_re_median
  }else{
    ResultTable_re_minmax <- merge(ResultTable_re_max, ResultTable_re_min, by = "Spp_PFAA")
    ResultTable_re <- merge(ResultTable_re_median, ResultTable_re_minmax, by = "Spp_PFAA")
  }
  
  return(ResultTable_re)
    
}


