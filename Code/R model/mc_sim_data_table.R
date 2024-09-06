

#' Title
#'
#' @param location 
#' @param season 
#' @param species_list 
#' @param PFAAs 
#' @param n_iter 
#' @param settings_obsDiet 
#' @param settings_modelDiet 
#' @param sediment_ratio_range 
#' @param sediment_ratio_only 
#'
#' @return
#' @export
#'
#' @examples
#' 
mcmc_input_files<-function(location,
                           season,
                           species_list,
                           PFAAs,
                           n_iter,
                           settings_obsDiet = NULL,
                           settings_modelDiet = NULL,
                           sediment_ratio_range = NULL,
                           sediment_ratio_only = FALSE){ # use a c(0,10) for a range of sediment ratios
  
  for (n_i in 1:n_iter){ # n_i is the number 'ith' in the iteration, random run process. 
  
  # 1. temperature -----
  temp<-rtnorm(n = 1000,
         mean = dwater.MEAN[c(dwater.MEAN$loc == location & 
                                              dwater.MEAN$Seasonality == season ), "Temp"],
         lower = as.numeric(dwater.MIN[c(dwater.MIN$loc == location &
                                              dwater.MIN$Seasonality == season), "Temp"]),
         upper = as.numeric(dwater.MAX[c(dwater.MAX$loc == location &
                                              dwater.MAX$Seasonality == season), "Temp"]),
         sd = (as.numeric(dwater.MIN[c(dwater.MIN$loc == location &
                                              dwater.MIN$Seasonality == season), "Temp"]) +
                 as.numeric(dwater.MAX[c(dwater.MAX$loc == location &
                                              dwater.MAX$Seasonality == season), "Temp"])) / 2)
  # hist(temp, breaks = 50, freq = TRUE, main = "temperature")
  # pick one random sample:
  temp_model<-sample(temp,1, replace=T) 
  if(sediment_ratio_only){
    temp_model<-as.numeric(dwater.MEAN[c(dwater.MEAN$loc == location & 
                                   dwater.MEAN$Seasonality == season ), "Temp"] )
  }
  
  # 2. species_list tissue concentrations -----
  for (i in 1:length(species_list)){
    for(j in 1:length(PFAAs)){
       pfaa<-rtnorm(n = 1000,
           mean = dfish.MEAN[c(dfish.MEAN$loc == location & 
                                                dfish.MEAN$Seasonality == season & 
                                                dfish.MEAN$SppAlias == species_list[i]), PFAAs[j]],
           lower = as.numeric(dfish.MIN[c(dfish.MIN$loc == location &
                                                dfish.MIN$Seasonality == season &
                                                dfish.MIN$SppAlias == species_list[i]), PFAAs[j]]),
           upper = as.numeric(dfish.MAX[c(dfish.MAX$loc == location &
                                                dfish.MAX$Seasonality == season &
                                                dfish.MAX$SppAlias == species_list[i]), PFAAs[j]]), 
           sd = abs(dfish.MEAN[c(dfish.MEAN$loc == location & 
                                                dfish.MEAN$Seasonality == season & 
                                                dfish.MEAN$SppAlias == species_list[i]), PFAAs[j]])+0.00001)
       
       if(j == 1){
         # random sampling:
         # hist(pfaa, breaks = 50, freq = TRUE, main = PFAAs[j])
         if(sediment_ratio_only){
           pfaa0<-as.numeric(dfish.MEAN[c(dfish.MEAN$loc == location & 
                                                dfish.MEAN$Seasonality == season & 
                                                dfish.MEAN$SppAlias == species_list[i]), PFAAs[j]])
         }else{
           pfaa0<-sample(pfaa, 1, replace = T)
         }
         # print(pfaa0)
       }else{
         if(sediment_ratio_only){
           pfaa0<-append(pfaa0, as.numeric(dfish.MEAN[c(dfish.MEAN$loc == location & 
                                                dfish.MEAN$Seasonality == season & 
                                                dfish.MEAN$SppAlias == species_list[i]), PFAAs[j]]))
         }else{
          pfaa0<-append(pfaa0, sample(pfaa, 1, replace = T))
         }
       }
    }
    assign(paste(species_list[i],"_pfaa", sep = ""), pfaa0) 
  }
  
  # 3. body size of species_list -----
  for (i in 1:length(species_list)){
    size<-rtnorm(n = 1000,
           mean = d.fish.w[c(d.fish.w$loc == location & 
                                d.fish.w$SppAlias == species_list[i]), "mean_mass"],
           # sd = 0.1,
           lower = as.numeric(d.fish.w[c(d.fish.w$loc == location &
                              d.fish.w$SppAlias == species_list[i]), "min_mass"]),
           upper = as.numeric(d.fish.w[c(d.fish.w$loc == location &
                              d.fish.w$SppAlias == species_list[i]), "max_mass"]),
           sd = c(as.numeric(d.fish.w[c(d.fish.w$loc == location &
                              d.fish.w$SppAlias == species_list[i]), "min_mass"])+
             as.numeric(d.fish.w[c(d.fish.w$loc == location &
                              d.fish.w$SppAlias == species_list[i]), "max_mass"])) /2)
    
     if(i == 1){
       # random sampling:
       # hist(size, breaks = 50, freq = TRUE, main = "size")
       if(sediment_ratio_only){
         size_model0<-as.numeric(d.fish.w[c(d.fish.w$loc == location & 
                                d.fish.w$SppAlias == species_list[i]), "mean_mass"])
       }else{
         size_model0<-sample(size, 1, replace = T)
       }
     }else{
       if(sediment_ratio_only){
         size_model0<-append(size_model0, as.numeric(d.fish.w[c(d.fish.w$loc == location & 
                                d.fish.w$SppAlias == species_list[i]), "mean_mass"]))
       }else{
         size_model0<-append(size_model0,sample(size, 1, replace = T))
       }
     }
     assign("Species_size", size_model0) 
  }
  
  # 4. dissolved oxygen (?) ------
  # 5. water PFAS concentrations ------
  for(j in 1:length(PFAAs)){
     pfaa<-rtnorm(n = 1000,
         mean = dwater.MEAN[c(dwater.MEAN$loc == location & 
                                              dwater.MEAN$Seasonality == season), PFAAs[j]],
         sd = dwater.MEAN[c(dwater.MEAN$loc == location & 
                                              dwater.MEAN$Seasonality == season), PFAAs[j]]+0.00001,
         lower = as.numeric(dwater.MIN[c(dwater.MIN$loc == location &
                                              dwater.MIN$Seasonality == season), PFAAs[j]]),
         
         upper = as.numeric(dwater.MAX[c(dwater.MAX$loc == location &
                                              dwater.MAX$Seasonality == season), PFAAs[j]]))
     if(j == 1){
       # random sampling:
       # hist(pfaa, breaks = 50, freq = TRUE, main = PFAAs[j])
       pfaa0<-sample(pfaa, 1, replace = T)
     }else{
       pfaa0<-append(pfaa0, sample(pfaa, 1, replace = T))
     }
     assign("Water_pfaa", pfaa0) 
     
  }
  
  # 6. sediment PFAS concentrations ------
  for (i in 1:length(species_list)){
    for(j in 1:length(PFAAs)){
      if(!is.null(sediment_ratio_range) & sediment_ratio_only){
       pfaa <- Water_pfaa[[j]] / (round(runif(1,sediment_ratio_range[1], sediment_ratio_range[2]), 3))
       pfaa <- pfaa/1000 # water is originally reported in ng/L in Brown et al., which would be ng/kg for solids. divide this by 1000 to get to ng/g
       # randomly selected one ratio 
       message("Estimate sediment:water proportions")
      }else{
        pfaa<-rtnorm(n = 1000,
           mean = dsed.MEAN[c(dsed.MEAN$loc == location & 
                                                dsed.MEAN$Seasonality == season), PFAAs[j]],
           sd = dsed.MEAN[c(dsed.MEAN$loc == location & 
                                                dsed.MEAN$Seasonality == season), PFAAs[j]]+0.00001,
           lower = as.numeric(dsed.MIN[c(dsed.MIN$loc == location &
                                                dsed.MIN$Seasonality == season ), PFAAs[j]]),
           upper = as.numeric(dsed.MAX[c(dsed.MAX$loc == location &
                                                dsed.MAX$Seasonality == season ), PFAAs[j]]))
       }
      if(j == 1){
         # random sampling:
         # hist(pfaa, breaks = 50, freq = TRUE, main = PFAAs[j])
         pfaa0<-sample(pfaa, 1, replace = T)
       }else{
         pfaa0<-append(pfaa0, sample(pfaa, 1, replace = T))
       }
    }
    assign("Sed_pfaa", pfaa0) 
  }
  
  # 7. make data tables for each location -----
  if(location == "JBA" & season == "Spring"){
  
    inputFiles_list<-create_data_tables(
      species = c("Phy", "Kil"	,"Chu",	"Dac", "Dar",	"Min",	"Mad",	"Pum",	"Swa"),
      group_species = c("plant", "fish", "fish", "fish", "fish", "fish", "fish", "fish", "fish"),
      WB_kg = c(NA,	Species_size),
      m_O = c(1, 1, 1, 1, 1, 1, 1, 1, 1),
      GRF = c(0.8, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150),
      P_B = c(0.5, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15), 
      diet = list( #
          Kil	= Kil_pfaa, # same
          Chu	= Chu_pfaa, # same
          Dac	= Dac_pfaa,
          Dar	= Dar_pfaa,
          Min	= Min_pfaa,
          Mad	= Mad_pfaa,
          Pum	= Pum_pfaa,
          Swa	= Swa_pfaa), # order for PFAA
    foodWeb = list(
        Phy = c(1,0,0,0,0,0,0,0,0,0),
        Kil = c(0.4,0.6,0,0,0,0,0,0,0,0),
        Swa = c(0.3,0.7,0,0,0,0,0,0,0,0),
        Dac = c(0.5,0.4,0,0,0,0,0,0,0,0),
        Min = c(0.4,0.6,0,0,0,0,0,0,0,0),
        Mad = c(0.3,0.4,0,0.1,0.1,0.1,0,0,0,0),
        Dar = c(0.4,0.5,0,0.1,0,0,0,0,0,0),
        Pum = c(0.4,0.3,0.1,0,0.1,0.1,0,0,0,0),
        Chu = c(0.5,0.2,0,0.1,0.1,0.1,0,0,0,0)
      ),
      C_WTO_ng_mL = Water_pfaa/1000, # order is for c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")
      C_s_ng_g = Sed_pfaa, # order of PFAA
      C_OX = 9.331,
      T = temp_model,
      log_Koc = c("Brown etal", 2.3317871	, 2.891004, 2.3032831,	3.0329491,	6,	5) # c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")
    ) 
  }
  
  if(location == "JBA" & season == "Summer"){
  
    inputFiles_list<-create_data_tables(
      species = c("Phy",	"Dac",	"Dar",	"Min",	"Fal",	"Bas",	"Mad",	"Pum",	"Swa"),
      group_species = c("plant", "fish", "fish", "fish", "fish", "fish", "fish", "fish", "fish"),
      WB_kg = c(NA,	Species_size),
      m_O = c(1, 1, 1, 1, 1, 1, 1, 1, 1),
      GRF = c(0.8, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150, 0.00150),
      P_B = c(0.5, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15), 
      diet = list( #
          Dac	= Dac_pfaa,
          Dar	= Dar_pfaa,
          Min	= Min_pfaa,
          Fal	= Fal_pfaa,
          Bas	= Bas_pfaa,
          Mad	= Mad_pfaa,
          Pum	= Pum_pfaa,
          Swa	= Swa_pfaa),
      foodWeb = list(
          Phy = c(0,0,0,0,0,0,0,0,0,0),
          Bas = c(0.8,0.2,0,0,0,0,0,0,0,0),
          Swa = c(0.3,0.6,0,0,0,0,0,0,0,0),
          Dac = c(0.5,0.4,0,0,0,0,0,0,0,0),
          Min = c(0.4,0.6,0,0,0,0,0,0,0,0),
          Fal = c(0.2,0.3,0,0.1,0.1,0.1,0,0,0,0),
          Mad = c(0.2,0.3,0,0.1,0.1,0.1,0,0,0,0),
          Dar = c(0.2,0.3,0,0.1,0,0,0,0,0,0),
          Pum = c(0.2,0.2,0.1,0,0.1,0.1,0,0,0,0)
          ),
      C_WTO_ng_mL = Water_pfaa/1000, # order is for c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")
      C_s_ng_g = Sed_pfaa, # order of PFAA
      C_OX = 6.449,
      log_Koc = c("Brown etal", 2.3317871	, 2.891004, 2.3032831,	3.0329491,	6,	5), # c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")
     T = temp_model
    ) 
  }
  
  if(location == "WG" & season == "Fall"){

    inputFiles_list<-create_data_tables(
      species = c("Phy", "Pry", "Bgl", "Bas"),
      group_species = c("plant", "fish", "fish", "fish"),
      WB_kg = c(NA,	Species_size),
      m_O = c(1, 1, 1, 1),
      GRF = c(0.8, 0.00150, 0.00150, 0.00150),
      P_B = c(0.5, 0.15, 0.15, 0.15), 
      diet = list( #
          Pry = Pry_pfaa,
          Bgl = Bgl_pfaa,
          Bas = Bas_pfaa),
      foodWeb = list(
          Phy = c(1,0,0,0,0),
          Pry = c(0.7,0.3,0,0,0),
          Bgl = c(0.1,0.3,0.6,0,0),
          Bas = c(0.1,0.1,0.7,0.1,0)
          ),
      C_WTO_ng_mL = Water_pfaa/1000, # order is for c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")
      C_s_ng_g = Sed_pfaa, # order of PFAA
      C_OX = 8.0,
      T = temp_model,
      log_Koc = c("Brown etal", 2.3317871	, 2.891004, 2.3032831,	3.0329491,	6,	5) # c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA")
    ) 
  }
  

  # 8. MODEL RUN -----
  if(!is.null(settings_modelDiet)){
      print(inputFiles_list)
      AllData.0<-runEcosystemModel(settings = settings_modelDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR',
                                    "Gill_uptake","Dietary_uptake", "Sediment_uptake",
                                    "Gill_up%", "Diet_up%", 'Sed_up%',
                                    "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = FALSE,
                  max_diet_WTO_C_s = FALSE)
  }
  
  if(!is.null(settings_obsDiet)){
      AllData.0<-runEcosystemModel(settings = settings_obsDiet,
                  inputFiles_list = inputFiles_list,
                  parameterList = c('C_B','C_WTO','C_s','k1','k2','ke','kg','kd',
                                    'Total Elimination', 'RMR',
                                    "Gill_uptake","Dietary_uptake", "Sediment_uptake",
                                    "Gill_up%", "Diet_up%", 'Sed_up%',
                                    "G_V", "G_D"),
                  PFAA_List = c('PFHxS','PFOS','PFOA','PFNA','PFDA','PFUA'),
                  dietData = inputFiles_list$dietData,
                  kRTable = kRTable,
                  min_diet_WTO_C_s = FALSE,
                  max_diet_WTO_C_s = FALSE)
  }

  AllData.0$iteration<-n_i
  # message(paste("Run: ", print(AllData.0$iteration[1])))
  
  if(n_i == 1){
      AllData <- AllData.0
    }else{
      AllData <- rbind(AllData, AllData.0)
    }
  }

  return(AllData)
}


