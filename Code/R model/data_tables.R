

# create and visualize bioaccumulation data tables, include metabolism

create_data_tables<-function(
    species, 
    group_species, # inv, fish, plant
    WB_kg, 
    foodWeb,
    diet,
    C_WTO_ng_mL,
    C_s_ng_g,
    C_OX,
    T, # order of PFAA
    min_diet = NULL,
    max_diet = NULL,
    MMR = NULL,
    RMR = NULL,
    pct_FeedRate = NULL,
    m_O = NULL,
    nu_NB = NULL,
    nu_LB = NULL,
    nu_PB = NULL,
    nu_OB = NULL,
    nu_WB = NULL,
    epsilon_N = NULL,
    epsilon_L = NULL,
    epsilon_P = NULL,
    epsilon_O = NULL,
    epsilon_W = NULL,
    A  = 0.00006,
    B  = 5.50000,
    GRF  = NULL,
    scaling_exp = 0.65,
    feeding_exp = 0.85,
    P_B = -9999, # whole fish protein content (%), only used for Liang et al models
    
    # chemicalData 
    C_WTO_max_ng_mL = NULL,
    C_WTO_min_ng_mL = NULL,
    C_s_max_ng_g = NULL,
    C_s_min_ng_g = NULL,
    # defaults
    Phi_exp_chem = c(0.9200000, 0.4000000, 0.80000000,  0.60000000, 0.350, 0.25),# order of PFAA
    
    # oceanData
    OCS = 0.03188235, # kg/L # organic carbon content in sediment
    delta_OCS = -99.0, # kg/L # 
    pH = 8.0, # pH 
    Phi_exp_ocean = 1.0,  # fraction of chemical in dossolved state
    C_SS = -99.0,  # kg/L # Concentration of suspended solids (kg/L)
    
    # chemicalParams
    PFAA = c("PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUA"), 
    chemID = c(NA, 6, 1, 8, 9, 10, 11),
    # first column is the source of the value
    Log_Kow = c("COSMOtherm", 5.200000, 6.4300000, 5.300000000, 5.920000,  6.500,  7.1500), # order of PFAA
    pKa     = c("Armitage",    0.000000, 0.0000000, 1.000000000, 1.000000,  1.000,  1.0000),# order of PFAA
    log_Kpw = c("Allendorf",   4.810000, 4.6700000, 4.200000000, 4.320000,  4.730,  4.6000),# order of PFAA
    log_Dmw = c("drodge",      3.820000, 4.8800000, 3.510000000, 4.040000,  4.630,  5.2200),# order of PFAA
    E_D     = c("Goeritz",     0.558000, 0.7210000, 0.138000000, 0.522000,  0.650,  0.7500),# order of PFAA
    E_W     = c("Martin10",    0.000790, 0.0676567, 0.000676567, 0.004851,  0.037,  0.1532),# order of PFAA
    log_Koc = c("Munoz",       2.331787, 2.8910036, 2.303283149, 3.032949,  4.580,  4.9900)# order of PFAA
  ){


  # create organism data **************************************
  if(is.null(m_O)){
    m_O <- vector()
    m_O <- append(c(rep(1, length(species))), m_O)
  }else{
    if(!c(length(m_O) == length(species))) {
      errorCondition("vector for m_O must be the same length as species")
    }
  }
  
  if(is.null(nu_NB)){
    nu_NB <- vector()
    for (i in 1:length(species)){
      if(group_species[i] == "inv" | group_species[i] == "plant"){
       nu_NB <- append(nu_NB, 0.02) 
      } 
      if(group_species[i] == "fish"){
       nu_NB <- append(nu_NB, 0.04)
      } 
    }
  }else{
    if(!c(length(nu_NB) == length(species))) {
      error("vector for nu_NB must be the same length as species")
    }
  }

  if(is.null(nu_LB)){
    nu_LB <- vector()
    nu_LB <- append(c(rep(0.01, length(species))), nu_LB)
  }else{
    if(!c(length(nu_LB) == length(species))) {
      error("vector for nu_LB must be the same length as organisms")
    }
  }

  if(is.null(nu_PB)){
    nu_PB <- vector()
    for (i in 1:length(species)){
      if(group_species[i] == "plant"){
       nu_PB <- append(nu_PB, 0) 
      } 
      if(group_species[i] == "inv"){
       nu_PB <- append(nu_PB, 0.001) 
      } 
      if(group_species[i] == "fish"){
       nu_PB <- append(nu_PB, 0.003)
      } 
    }
  }else{
    if(!c(length(nu_PB) == length(species))) {
      error("vector for nu_PB must be the same length as species")
    }
  }
  
  if(is.null(nu_OB)){
    nu_OB <- vector()
    nu_OB <- append(c(rep(0.15, length(species))), nu_OB)
  }else{
    if(!c(length(nu_OB) == length(species))) {
      error("vector for nu_OB must be the same length as organisms")
    }
  }
  
  if(is.null(nu_WB)){
    nu_WB <- vector()
    for (i in 1:length(species)){
      if(group_species[i] == "inv" | group_species[i] == "plant"){
       nu_WB <- append(nu_WB, 0.819) 
      } 
      if(group_species[i] == "fish"){
       nu_WB <- append(nu_WB, 0.797)
      } 
    }
  }else{
    if(!c(length(nu_WB) == length(species))) {
      error("vector for nu_WB must be the same length as species")
    }
  }
  
  if(is.null(epsilon_N)){
    epsilon_N <- vector()
    for (i in 1:length(species)){
      if(group_species[i] == "plant"){
       epsilon_N <- append(epsilon_N, NA) 
      } 
      if(group_species[i] == "inv"){
       epsilon_N <- append(epsilon_N, 0.75) 
      } 
      if(group_species[i] == "fish"){
       epsilon_N <- append(epsilon_N, 0.92)
      } 
    }
  }else{
    if(!c(length(epsilon_N) == length(species))) {
      error("vector for epsilon_N must be the same length as species")
    }
  }

  if(is.null(epsilon_L)){
    epsilon_L <- vector()
    for (i in 1:length(species)){
      if(group_species[i] == "plant"){
       epsilon_L <- append(epsilon_L, NA) 
      } 
      if(group_species[i] == "inv"){
       epsilon_L <- append(epsilon_L, 0.75) 
      } 
      if(group_species[i] == "fish"){
       epsilon_L <- append(epsilon_L, 0.92)
      } 
    }
  }else{
    if(!c(length(epsilon_L) == length(species))) {
      error("vector for epsilon_L must be the same length as species")
    }
  }
  
  if(is.null(epsilon_P)){
    epsilon_P <- vector()
    for (i in 1:length(species)){
      if(group_species[i] == "plant"){
       epsilon_P <- append(epsilon_P, NA) 
      } 
      if(group_species[i] == "inv"){
       epsilon_P <- append(epsilon_P, 0.75) 
      } 
      if(group_species[i] == "fish"){
       epsilon_P <- append(epsilon_P, 0.92)
      } 
    }
  }else{
    if(!c(length(epsilon_P) == length(species))) {
      error("vector for epsilon_P must be the same length as species")
    }
  }
  
  if(is.null(epsilon_O)){
    epsilon_O <- vector()
    for (i in 1:length(species)){
      if(group_species[i] == "plant"){
       epsilon_O <- append(epsilon_O, NA) 
      } 
      if(group_species[i] == "fish" | group_species[i] == "inv"){
       epsilon_O <- append(epsilon_O, 0.6)
      } 
    }
  }else{
    if(!c(length(epsilon_O) == length(species))) {
      error("vector for epsilon_O must be the same length as species")
    }
  }

  if(is.null(epsilon_W)){
    epsilon_W <- vector()
    for (i in 1:length(species)){
      if(group_species[i] == "plant"){
       epsilon_W <- append(epsilon_W, NA) 
      } 
      if(group_species[i] == "fish" | group_species[i] == "inv"){
       epsilon_W <- append(epsilon_W, 0.7)
      } 
    }
  }else{
    if(!c(length(epsilon_W) == length(species))) {
      error("vector for epsilon_W must be the same length as species")
    }
  }
  
  if(is.null(A)){
    A <- c(0.00006, rep(NA, length(species)-1)) 
  }else{
    A <- c(A, rep(NA, length(species)-1))
  }

  if(is.null(B)){
    B <- c(5.5, rep(NA, length(species)-1)) 
  }else{
    B <- c(B, rep(NA, length(species)-1))
  }
  
  if(is.null(GRF)){
    GRF <- vector()
    for (i in 1:length(species)){
      if(group_species[i] == "plant"){
       GRF <- append(GRF, 0.08) 
      } 
      if(group_species[i] == "inv"){
       GRF <- append(GRF, 0.00035)
      } 
      if(group_species[i] == "fish"){
       GRF <- append(GRF, 0.0014)
      } 
    }
  }else{
    if(!c(length(GRF) == length(species))) {
      error("vector for GRF must be the same length as species")
    }
  }
  
  if(is.null(MMR)){
    MMR <- vector()
    for (i in 1:length(species)){
      if(group_species[i] == "plant"){
       MMR <- append(MMR, NA) 
      } 
      if(group_species[i] == "inv"){
       MMR <- append(MMR, NA)
      } 
      if(group_species[i] == "fish"){
       MMR <- append(MMR, NA)
      } 
    }
  }else{
    if(!c(length(MMR) == length(species))) {
      error("vector for MMR must be the same length as species")
    }
  }
  if(is.null(RMR)){
    RMR <- vector()
    for (i in 1:length(species)){
      if(group_species[i] == "plant"){
       RMR <- append(RMR, NA) 
      } 
      if(group_species[i] == "inv"){
       RMR <- append(RMR, NA)
      } 
      if(group_species[i] == "fish"){
       RMR <- append(RMR, NA)
      } 
    }
  }else{
    if(!c(length(RMR) == length(species))) {
      error("vector for RMR must be the same length as species")
    }
  }
  if(is.null(pct_FeedRate)){
    pct_FeedRate <- vector()
    for (i in 1:length(species)){
      if(group_species[i] == "plant"){
       pct_FeedRate <- append(pct_FeedRate, NA) 
      } 
      if(group_species[i] == "inv"){
       pct_FeedRate <- append(pct_FeedRate, -999)
      } 
      if(group_species[i] == "fish"){
       pct_FeedRate <- append(pct_FeedRate, 4)
      } 
    }
  }else{
    if(!c(length(pct_FeedRate) == length(species))) {
      error("vector for pct_FeedRate must be the same length as species")
    }
  }
  
  switchk_1<-vector()
  for (i in 1:length(species)){
    if(group_species[i] == "plant"){
     switchk_1 <- append(switchk_1, 0) 
    } 
    if(group_species[i] == "inv"){
     switchk_1 <- append(switchk_1, 1)
    } 
    if(group_species[i] == "fish"){
     switchk_1 <- append(switchk_1, 1)
    } 
  }
  switchG_D<-vector()
  for (i in 1:length(species)){
    if(group_species[i] == "plant"){
     switchG_D <- append(switchG_D, 0) 
    } 
    if(group_species[i] == "inv"){
     switchG_D <- append(switchG_D, 2)
    } 
    if(group_species[i] == "fish"){
     switchG_D <- append(switchG_D, 1)
    } 
  }
  switchk_R<-vector()
  for (i in 1:length(species)){
    if(group_species[i] == "plant"){
     switchk_R <- append(switchk_R, 0) 
    } 
    if(group_species[i] == "inv"){
     switchk_R <- append(switchk_R, 0)
    } 
    if(group_species[i] == "fish"){
     switchk_R <- append(switchk_R, 1)
    } 
  }
  
  if(length(scaling_exp) == 1){
    scaling_exp <- rep(scaling_exp, length(species))
  }
  if(length(feeding_exp) == 1){
    feeding_exp <- rep(feeding_exp, length(species))
  }
  if(length(P_B) == 1){
    P_B <- rep(P_B, length(species))
  }else{
    P_B <- P_B # value for each species
  }

  
  organismData <- as.data.frame(matrix(data = 
                                         c(WB_kg, m_O, nu_NB, nu_LB, nu_PB, nu_OB,
                                           nu_WB, epsilon_N, epsilon_L, epsilon_P,
                                           epsilon_O, epsilon_W, A, B,
                                           GRF, switchk_1, switchG_D, switchk_R,
                                           MMR, RMR, pct_FeedRate, scaling_exp,
                                           feeding_exp, P_B),
                                       byrow = TRUE, ncol = length(species)))

  colnames(organismData) <- species
  rownames(organismData) <- c("W_B", "m_O", "nu_NB", "nu_LB", "nu_PB", "nu_OB", 
                            "nu_WB", "epsilon_N", "epsilon_L", "epsilon_P",
                            "epsilon_O", "epsilon_W", "A", "B",
                            "GRF", "switchk_1", "switchG_D", "switchk_R",
                            "MMR", "RMR", "pct_FeedRate", "scaling_exp", "feeding_exp", "P_B")


  # 2. chemicalData **********************
  # chemicalData 
  # order of PFAA
  if(is.null(C_WTO_max_ng_mL)){
    C_WTO_max_ng_mL <- C_WTO_ng_mL * 1.1 # 10 % higher 
  }
  if(is.null(C_WTO_min_ng_mL)){
    C_WTO_min_ng_mL <- C_WTO_ng_mL * 0.9 # 10 % lower
    message("min and max PFAS in water: 10 % +/- water PFAS")
  }
  if(is.null(C_s_max_ng_g)){
    C_s_max_ng_g <- C_s_ng_g * 1.1 # 10 % higher 
  }
  if(is.null(C_s_min_ng_g)){
    C_s_min_ng_g <- C_s_ng_g * 0.9 # 10 % lower
    message("min and max PFAS in sediment: 10 % +/- sediment PFAS")
  }
  

  
  chemicalData <- data.frame(matrix(data = c(C_WTO_ng_mL, C_WTO_max_ng_mL, C_WTO_min_ng_mL,
                                             C_s_ng_g, C_s_max_ng_g, C_s_min_ng_g, Phi_exp_chem),
                                    ncol=length(PFAA), byrow=TRUE))
  colnames(chemicalData) <-c(PFAA)
  rownames(chemicalData) <-c("C_WTO", "C_WTO_max", "C_WTO_min", 
                             "C_s", "C_s_max", "C_s_min", "Phi_exp")
  chemicalData$units_def<-c("ng/mL",
                            "ng/mL",
                            "ng/mL",
                            "ng/g (dw)",
                            "ng/g (dw)",
                            "ng/g (dw)",
                            "unitless")
  
  # 3. oceanData *****************   # units need?  here here 
  oceanData <- data.frame(matrix(data = c(C_OX, T, OCS, delta_OCS, pH, Phi_exp_ocean, C_SS),
                                 ncol=length(1), byrow=TRUE))
  colnames(oceanData) <-c("value")
  rownames(oceanData) <-c("C_OX", "T", "OCS", "delta_OCS", "pH", "Phi_exp", "C_SS")
  oceanData$units_def<-c("mg/L",
                         "ºC",
                         "kg/L",
                         "kg/L",
                         "unitless",
                         "unitless",
                         "n/a?")

  # 4. chemicalParams *****************
  # the order of PFAA
  chemicalParams <- data.frame(matrix(data = c(Log_Kow, pKa, log_Kpw, log_Dmw,
                                               E_D, E_W, log_Koc, chemID),
                                 ncol=length(PFAA)+1, byrow=TRUE))
  colnames(chemicalParams) <-c("sourceVal", PFAA)
  rownames(chemicalParams) <-c("Log_Kow", "pKa", "log_Kpw", "log_Dmw", "E_D", "E_W", "log_Koc", "chemID")
  chemicalParams$units<-c("unitless", 
                          "unitless",
                          "unitless",
                          "unitless",
                          "unitless", 
                          "unitless",
                          "Koc = Kd / OCS",
                          "n")
  
  
  # 5. foodWebData *****************
  # how much each species eat what
  # fill out each row, in a list 
  for (i in 1:length(foodWeb)){
    if( i == 1){
      foodWebData<-foodWeb[[1]]
    }else{
      foodWebData<-rbind(foodWebData, foodWeb[[i]])
    }
  }
  rownames(foodWebData)<-names(foodWeb)
  colnames(foodWebData)<-c("Sed", names( foodWeb))
  
  # 6. dietData median, min, max diet. 
  # the order of PFAA
  for (i in 1:length(diet)){
    if( i == 1){
      dietData<-diet[[1]]
    }else{
      dietData<-rbind(dietData, diet[[i]])
    }
  }
  rownames(dietData)<-names(diet)
  colnames(dietData)<-PFAA
  
  if(!is.null(min_diet)){
    for (i in 1:length(min_diet)){
      if( i == 1){
        min_dietData<-min_diet[[1]]
      }else{
        min_dietData<-rbind(min_dietData, min_diet[[i]])
      }
    }
    rownames(min_dietData)<-names(min_diet)
    colnames(min_dietData)<-PFAA
    min_dietData<-as.data.frame(min_dietData)
  }else{
    min_diet = "Not available"
  }

  
  if(!is.null(max_diet)){
    for (i in 1:length(max_diet)){
      if( i == 1){
        max_dietData<-max_diet[[1]]
      }else{
        max_dietData<-rbind(max_dietData, max_diet[[i]])
      }
    }
    rownames(max_dietData)<-names(max_diet)
    colnames(max_dietData)<-PFAA
    max_dietData<-as.data.frame(max_dietData)
  }else{
    max_diet = "Not available"
  }

  
  chemicalParams$PFHxS<-as.numeric(chemicalParams$PFHxS)
  chemicalParams$PFOS<-as.numeric(chemicalParams$PFOS)
  chemicalParams$PFOA<-as.numeric(chemicalParams$PFOA)
  chemicalParams$PFNA<-as.numeric(chemicalParams$PFNA)
  chemicalParams$PFDA<-as.numeric(chemicalParams$PFDA)
  chemicalParams$PFUA<-as.numeric(chemicalParams$PFUA)
  
  numSpecies<-length(species)# assuming sediment is always there 
  foodWebData<-as.data.frame(foodWebData)
  dietData<-as.data.frame(dietData)
  
  if(length(min_diet) == 1 & length(max_diet) == 1){
    inputFiles_list = tibble::lst(numSpecies, 
                          organismData,
                           chemicalData,
                           oceanData,
                           chemicalParams,
                           foodWebData,
                           dietData)

  } else {
    inputFiles_list = tibble::lst(numSpecies, 
                              organismData,
                               chemicalData,
                               oceanData,
                               chemicalParams,
                               foodWebData,
                               dietData,
                               min_dietData,
                               max_dietData)
  }

  return(inputFiles_list)
}

  
  
  
  
  
  
