
#' Title: Food Web Model - run for one PFAA at the time. 
#'
#' @param PFAA 
#' @param settings 
#' @param numSpecies 
#' @param oceanData 
#' @param chemicalData 
#' @param chemicalParams 
#' @param organismData 
#' @param foodWebData 
#' @param dietData 
#' @param kRTable 
#'
#' @return
#' @export
#'
#' @examples
BioaccumulationModel <- function(
    PFAA, # a string c()
    settings, # a list new_Settings()
    numSpecies,# numeric, integer
    oceanData, # DataFrame  
    chemicalData, # DataFrame 
    chemicalParams, # DataFrame 
    organismData, # DataFrame 
    foodWebData, # pd.DataFrame 
    dietData = NULL, # Optional[pd.DataFrame]=None,
    kRTable = NULL # :Optional[pd.DataFrame]=None
  ){

    ###################################*
    # Set up model settings ------
    ###################################*
    chemicalData = chemicalData[,PFAA, drop = FALSE] # column with select PFAA, all rows

    chemicalParams = chemicalParams[, PFAA, drop = FALSE] # column with select PFAA, all rows
    
    if (!is.null(dietData)){
        dietData = dietData[, PFAA, drop = FALSE] # column with select PFAA, all rows
    }

    #############################*
    # Read in parameters -------
    #############################*

    ## Ocean Parameters (numeric; feed into new_Environment())
    C_OX = oceanData['C_OX',1] #Dissolved Oxygen Concentration (mg/L)
    T = oceanData['T',1] #Temperature (C)
    C_SS = oceanData['C_SS', 1] #Concentration of suspended solids (kg/L)
    OCS = oceanData['OCS', 1] #Organic carbon content of the sediment (fraction)
    pH = oceanData['pH',1] # ocean pH
    pHi = 7.4 # default for internal pH as suggested by Armitage et  al. 2013. # also in class 
    pHg = pH - 1 # also in class 
    
    ## Chemical Parameters (numeric; feed into new_Chemical)
    Log_Kow = chemicalParams['Log_Kow', ] # the whole row, all PFAA
    Log_Kpw = log10(1.36 * 10^chemicalParams[ 'log_Kpw', ])
    Log_Dmw = chemicalParams['log_Dmw', ] # Dmw measured by Droge 2019
    pKa = chemicalParams[ 'pKa', ]
    Log_Koc = chemicalParams['log_Koc', ] # choose values relevant for model system
    E_W_exp = chemicalParams['E_W', ] # calculated at Cox = 10 mg/L. Estimated based on T=12, DO sat approximately 92%
    E_D_exp = chemicalParams['E_D', ] # new for v2

    # chemID
    chemID = chemicalParams['chemID', ]
    # print(chemID)
    
    ## Organism Parameters
    W_B = organismData['W_B',] #Final Weight of the organism (kg)
    m_O = organismData['m_O',] #Fraction of the respiratory ventilation that involves overlying water

    nu_NB = organismData['nu_NB',]  #Neutral lipid fraction in the gut 
    nu_LB = organismData['nu_LB',]  #Phosphoipid fraction in the gut 
    nu_PB = organismData['nu_PB',]  #Protein fraction in the gut 
    nu_OB = organismData['nu_OB',]  #NLOM (or "other") fraction. in the gut 
    nu_WB = organismData['nu_WB',]  #Water content in the gut 
    
    P_B = organismData['P_B',] 

    epsilon_N = organismData['epsilon_N',]  #Dietary assimilation efficiency of neutral lipid contents of diet (fraction)
    epsilon_L = organismData['epsilon_L',]  #Dietary assimilation efficiency of phospholipid contents of diet ((fraction)
    epsilon_P = organismData['epsilon_P',]  #Dietary assimilation efficiency of protein contents of diet (fraction)
    epsilon_O = organismData['epsilon_O',]  #Dietary assimilation efficiency of the NLOM content of diet (fraction)
    epsilon_W = organismData['epsilon_W',]  #Dietary assimilation efficiency of water contents of diet (fraction)

    A = organismData['A',]  #Aqueous phase resistance constant for k1 calculation
    B = organismData['B',]  #Organic phase resistance constant for k1 calculation
    GRF = organismData['GRF',] #growth weight factor (for growth dilution calculations)

    switchk_1 = organismData['switchk_1',] #dependent on organism type
    switchG_D = organismData['switchG_D',]  #dependent on organism type
    switchk_R = organismData['switchk_R',] # dependent on organism type
    
    MMR = organismData['MMR',]
    RMR = organismData['RMR',]
    scaling_exp = organismData['scaling_exp',]
    feeding_exp = organismData['feeding_exp',]
    
    sigma =  1 # sigma is the scavenging efficiency of particles
    Beta =  0.035 # Proportionality constant representing the sorption capacity of NLOM to that of octanol.
    # Beta is specifically relevant to calculating partitioning in phytoplankton, where NLOM is replaced with NLOC

    #############################################*
    # Read in system-specific parameters --------
    #############################################*
    ## System-Specific Chemical Data
    C_WTO = chemicalData['C_WTO', ] #Total chemical concentration in the water column above the sediments (g/L)
    C_s = chemicalData['C_s', ] #.0358; %C_s is the Concentration in sediment
    Phi = chemicalData['Phi_exp', ] # Field-measured Phi reported for a given system

    # Food Web Data
    Pd = foodWebData[,1] #Fraction of the diet consisting of detritus
    P = foodWebData[,c(1:numSpecies+1), drop = FALSE] #Fraction of the diet consisting of prey item i

    ########################################*
    # Calculate dietary composition ------
    ########################################*

    nu_ND = array(rep(0, numSpecies))
    nu_LD = array(rep(0, numSpecies))
    nu_PD = array(rep(0, numSpecies))
    nu_OD = array(rep(0, numSpecies))
    nu_WD = array(rep(0, numSpecies))

    for (i in 1:numSpecies){
          nu_ND[i] = sum((P[i,] * nu_NB)) #Overall neutral lipid content of diet
          nu_LD[i] = sum((P[i,] * nu_LB)) #Overall phospholipid content of the diet
          nu_PD[i] = sum((P[i,] * nu_PB)) #Overall protein content of the diet
          nu_OD[i] = sum((P[i,] * nu_OB)) #Overall NLOM content of the diet
          nu_WD[i] = sum((P[i,] * nu_WB)) #Overall water content of diet
    } 
     
    ###############################*
    # Instantiate classes --------
    ###############################*
    env = new_Environment(C_OX=C_OX, T=T, C_SS=C_SS, OCS=OCS, pH=pH)$get_vars()
    chem = new_Chemical(Log_Kow=Log_Kow, Log_Kpw=Log_Kpw, Log_Dmw=Log_Dmw, pKa=pKa,
                        Log_Koc=Log_Koc, chemID=chemID,
                        Ed_exp=E_D_exp, Ew_exp=E_W_exp, 
                        en =env,
                        settings=settings)$get_vars()
    
    chemdata = new_ChemData(C_WTO=C_WTO, C_s=C_s, Phi=Phi, env=env, chem=chem)$get_vars()

    #####################################*
    # Run model and save results --------
    #####################################*
  
    Parameters = array(data = c('C_B','Phi','mO','C_WTO','C_WDP',
                                'Water','C_s','FeedRate', 'Sediment', 'Diet', 'Gill_uptake',
                                'Dietary_uptake', 'Sediment_uptake','Gill_up%','Diet_up%', 'Sed_up%','G_V','G_D',
                                'G_F','W_B','Ew','Ed', 'k1',
                                'k2', 'kd','ke','kg','kr_est',
                                'Total Elimination','kr_pct','pKa','logDbw','log Kow',
                                'log Dmw','log Dow','log Kpw','D_BW','D_MW',
                                'D_OW','K_PW','pHi','pHg','K_GB',
                                'nu_NB','nu_LB','nu_PB','nu_OB','nu_WB',
                                'epsilon_N','epsilon_L','epsilon_P', 'epsilon_O','Log_Koc',
                                'Phi','RMR'))
    
    Results = as.data.frame(matrix(0, numSpecies, length(Parameters)))
    colnames(Results)<-Parameters
    rownames(Results)<-colnames(organismData)
  
    if(settings$chooseDiet == 'default'){
        C_D = array(rep(0, numSpecies))

        for (i in 1:numSpecies){ # not phytoplankton, start with numSpecies 1: 4
            orgi = new_Organism(W_B = W_B[i], m_O = m_O[i], GRF = GRF[i], 
                                nu_NB = nu_NB[i], nu_LB = nu_LB[i],
                                nu_PB = nu_PB[i], nu_OB = nu_OB[i],
                                nu_WB = nu_WB[i], nu_ND = nu_ND[i],
                                nu_LD = nu_LD[i], nu_PD = nu_PD[i],
                                nu_OD =nu_OD[i], nu_WD = nu_WD[i],
                                epsilon_N = epsilon_N[i], epsilon_L = epsilon_L[i],
                                epsilon_P = epsilon_P[i], epsilon_O = epsilon_O[i],
                                epsilon_W = epsilon_W[i], switchk_1 = switchk_1[i],
                                switchG_D = switchG_D[i], switchk_R = switchk_R[i],
                                A = A[i], B = B[i],
                                MMR = MMR[i], RMR = RMR[i],
                                P_B = P_B[i],
                                scaling_exp = scaling_exp[i],
                                feeding_exp = feeding_exp[i],
                                sigma = sigma,
                                env = env, settings = settings,chem = chem)$get_vars()
            Results0<-SSC_B(settings=settings, chemdata=chemdata, C_D=C_D, P=P[i,], Pd=Pd[i], env=env, chem=chem, org=orgi, kRTable = kRTable)
            Results[i,] = Results0
            C_D[i] = Results[i,1]
            # print(C_D)
            # print("**")
        }

    } else if (settings$chooseDiet == 'forced'){
        # read in empirical diet data
        C_D <- dietData # manually provided PFAS in each diet item
        
        # insert modeled phytoplankton data, which are not empirically measured
        org0 = new_Organism(W_B[1], m_O[1], GRF[1], nu_NB[1], nu_LB[1], nu_PB[1], nu_OB[1], nu_WB[1],
                        nu_ND[1], nu_LD[1], nu_PD[1], nu_OD[1], nu_WD[1], epsilon_N[1], epsilon_L[1], epsilon_P[1],
                        epsilon_O[1], epsilon_W[1], switchk_1[1], switchG_D[1], switchk_R[1], A[1], B[1], sigma,
                        MMR = MMR[1], RMR = RMR[1], P_B = P_B[1],
                        scaling_exp = scaling_exp[1],feeding_exp = feeding_exp[1],
                        env = env, settings = settings,chem = chem)$get_vars()
        # Results.iloc[0,:] = getSSC_B(settings, chemdata, C_D, P.iloc[i,:], Pd[i], env, chem, org0,kRTable)
        
        Results_phy<-SSC_B(settings=settings, chemdata=chemdata,
                           C_D=C_D, P=P[1,], Pd=Pd[1], env=env, chem=chem, org=org0, kRTable = kRTable)

        Results[1,] = Results_phy
        C_D[1,1] = Results_phy[1,1] # new estimated value, used in the loop below
      
        # print(numSpecies)
        for (i in 1:numSpecies){ # all including phytoplankton
            orgi = new_Organism(W_B[i], m_O[i], GRF[i], nu_NB[i], nu_LB[i], nu_PB[i], nu_OB[i], nu_WB[i],
                            nu_ND[i], nu_LD[i], nu_PD[i], nu_OD[i], nu_WD[i], epsilon_N[i], epsilon_L[i], epsilon_P[i],
                            epsilon_O[i], epsilon_W[i], switchk_1[i], switchG_D[i], switchk_R[i], A[i], B[i], sigma,
                            MMR = MMR[i], RMR = RMR[i], P_B = P_B[i],
                            scaling_exp = scaling_exp[i],feeding_exp = feeding_exp[i],
                            env = env, settings = settings,chem = chem)$get_vars()
            Results0 = SSC_B(settings=settings, chemdata=chemdata,
                             C_D=C_D, P=P[i,], Pd=Pd[i], env=env, chem=chem, org=orgi, kRTable = kRTable)
            Results[i,] = Results0
        }
    }else if (settings$chooseDiet == 'zero'){
      
       C_D = array(rep(0, numSpecies))

        for (i in 1:numSpecies){ # not phytoplankton, start with numSpecies 1: 4
            orgi = new_Organism(W_B = W_B[i], m_O = m_O[i], GRF = GRF[i], 
                                nu_NB = nu_NB[i], nu_LB = nu_LB[i],
                                nu_PB = nu_PB[i], nu_OB = nu_OB[i],
                                nu_WB = nu_WB[i], nu_ND = nu_ND[i],
                                nu_LD = nu_LD[i], nu_PD = nu_PD[i],
                                nu_OD =nu_OD[i], nu_WD = nu_WD[i],
                                epsilon_N = epsilon_N[i], epsilon_L = epsilon_L[i],
                                epsilon_P = epsilon_P[i], epsilon_O = epsilon_O[i],
                                epsilon_W = epsilon_W[i], switchk_1 = switchk_1[i],
                                switchG_D = switchG_D[i], switchk_R = switchk_R[i],
                                A = A[i], B = B[i],
                                MMR = MMR[i], RMR = RMR[i],
                                P_B = P_B[i],
                                scaling_exp = scaling_exp[i],
                                feeding_exp = feeding_exp[i],
                                sigma = sigma,
                                env = env, settings = settings,chem = chem)$get_vars()
            Results0<-SSC_B(settings=settings, chemdata=chemdata, C_D=C_D, P=P[i,], Pd=Pd[i], env=env, chem=chem, org=orgi, kRTable = kRTable)
            Results[i,] = Results0
            # C_D[i] = Results[i,1]
            # print(C_D)
            # print("**")
        }
    }
    
    Results$C_B <- Results$C_B*1000 # convert from ng/g to ng/kg
    Results$logC_B <- log10(Results$C_B) # units in ng/kg
    Results$C_B_ngg <- Results$C_B / 1000 # convert back to ng/g
    Results$logBAF <- log10(Results$C_B_ngg / Results$C_WTO)
    Results$BMF <- Results$C_B_ngg / Results$Diet
    
    return(Results)

}

     



  
    
