  
#' Title: Steady State Equation for calculating concentration in a single organism
#'
#' @param settings 
#' @param PFAA 
#' @param chemdata 
#' @param env 
#' @param chem 
#' @param org 
#' @param C_D 
#' @param P 
#' @param Pd 
#' @param kRTable 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
SSC_B<-function(settings,
                PFAA,
                chemdata,
                env,
                chem,
                org, 
                C_D,
                P,
                Pd,
                kRTable=NULL){
  
    ########################################################### *
    ## Calculate intermediate model inputs & rates ----
    ########################################################### *
    C_WDP <- chemdata$C_WDP

    k_1 <- org$k_1 #k_1 is the clearance rate constant (L/kg*d) for chemical uptake via respiratory area
    
    if (settings$chooseModel == "Sun_etal_2022"){
      k_2 = k_1 / org$D_BW  
      k_M = 0
    }else if(settings$chooseModel == "Liang_etal_2022"){
      k_2 = k_1 / (org$P_B * chem$K_OW)    
      k_M = 0 
  
      # k_BTR = k_2 + k_E + k_G # reference non biotranformed chemical 
      # k_BT = k_2 + k_E + k_G + k_M # biotransformed chemical 
      # k_M = k_BT - k_BTR # new 
      
    }else{
      message("expect error b/c the model type is non specified")
    }
    
    if(org$switchk_1 == 0){ # Phytoplankton
        k_D = 0
        k_E = 0
        k_G = org$GRF
    } else if (org$switchk_1 == 1){ # Zooplankton, Aquatic Invertebrates, Fish
        k_D = org$k_D
        k_E = org$k_E  #k_E is the rate constant (1/d) for chemical elimination via excretion into egested feces
        k_G = org$k_G  #k_G is the growth rate constant
    } else {
        stop('error in k1 switch selection for dietary & growth uptake & elimination pathways')
    }

    if (org$switchk_R == 1) {
        # df_krkbRatio
        k_R_est = kRTable[kRTable$chemID == chem$chemID,'kr/kb'] * k_2
    } else if (org$switchk_R == 0) {
        k_R_est = 0
    } else { 
        stop('error in kR switch')
    }
    
    
     ############################################### *
     ## Calculate tissue concentration -----
     ############################################### *
    
    if (settings$chooseModel == "Sun_etal_2022"){
      C_B = (( (k_1 * (org$m_O * chemdata$Phi * chemdata$C_WTO + (1 - org$m_O) * C_WDP))  +
              (k_D * (sum(P * as.numeric(unlist(C_D))) + (Pd * chemdata$C_s))) ) /
              (k_2 + k_E + k_M + k_G + k_R_est))

    }else if(settings$chooseModel == "Liang_etal_2022"){
      C_B = (k_1 * chemdata$C_WTO) + 
            (k_D * (sum(P * as.numeric(unlist(C_D)))))
      C_B = C_B - (k_2 + k_E + k_M + k_G) * C_B
      
    }else{
      message("expect error b/c model type")
    }
    
    ################################################ *
    ## Create table of output values ------
    ################################################ *
    D_OW = chem$D_OW
    D_MW = chem$D_MW
    D_BW = org$D_BW

    G_F = org$G_F
    #G_F = 0
    G_D = org$G_D
    #G_D = 1
    G_V = org$G_V
    E_D = chem$E_D
    K_GB = org$K_GB
    E_W = chem$E_W

    # *************** Organism specific chemical uptake here ***********
    # Diet = sum(P * t(C_D)) + (Pd * chemdata$C_s) # full diet 
    # split up sediment and diet 
    Diet = sum(P * t(C_D)) # only food 
    Sediment = (Pd * chemdata$C_s) # only sediment
    Water = (org$m_O * chemdata$Phi * chemdata$C_WTO + (1-org$m_O) * C_WDP)
    
    FeedRate = G_D / org$W_B
    
              # ((k_1 * (org$m_O * chemdata$Phi * chemdata$C_WTO + (1 - org$m_O) * C_WDP) 
    Gill_uptake = k_1 * (org$m_O * chemdata$Phi * chemdata$C_WTO + (1 - org$m_O) * C_WDP) # L/kg*d * ng/mL (or g/L) = g chemical/kg fish/day
    # Dietary_uptake = k_D * (sum(P * t(C_D)) + (Pd * chemdata$C_s)) # kg food/kg org * ng/g (or g/kg food) = g chemical/kg fish/day
    
    # split up sediment and diet 
                   # k_D * (sum(P * as.numeric(unlist(C_D))) + (Pd * chemdata$C_s)))
    Dietary_uptake = k_D * (sum(P * t(C_D))) # kg food/kg org * ng/g (or g/kg food) = g chemical/kg fish/day
    Sediment_uptake = k_D * (Pd * chemdata$C_s) # kg food/kg org * ng/g (or g/kg food) = g chemical/kg fish/day

    Uptake = Gill_uptake + Dietary_uptake + Sediment_uptake
    Gill_uppct = Gill_uptake / Uptake
    Diet_uppct = Dietary_uptake / Uptake
    Sediment_uppct = Sediment_uptake / Uptake
    TotalElim_rate = (k_2 + k_E + k_M + k_G + k_R_est)  
    # *************** Organism specific chemical uptake ***********

    NL_pct = (org$nu_NB * D_OW) / D_BW
    PL_pct = (org$nu_LB * D_MW) / D_BW
    Protein_pct = (org$nu_PB * chem$K_PW) / D_BW # protein here 
    NLOM_pct = (org$nu_OB * D_OW * 0.05) / D_BW
    water_pct = (org$nu_WB) / D_BW

    kr_pct = k_R_est / (k_2 + k_E + k_G + k_R_est)

    Output_Data = t(as.data.frame(c("C_B" = C_B,
                                  'mO' = org$m_O, 
                                  'Phi' = chemdata$Phi,
                                  'C_WTO' = chemdata$C_WTO,
                                  'C_WDP' = C_WDP,
                                  'Water' = Water,
                                  'C_s' = chemdata$C_s,
                                  'FeedRate' = FeedRate,
                                  'Sediment' = Sediment,
                                  'Diet' = Diet,
                                  'Gill_uptake' = Gill_uptake,
                                  'Dietary_uptake' = Dietary_uptake,
                                  'Sediment_uptake' = Sediment_uptake,
                                  'Gill_up' = Gill_uppct,
                                  'Diet_up' = Diet_uppct,
                                  'Sediment_up' = Sediment_uppct,
                                  'G_V' = G_V,
                                  'G_D' = G_D,
                            'G_F' = G_F,
                            'W_B' = org$W_B,
                            'Ew' =  E_W,
                            'Ed' = E_D,
                            'k1' = k_1,
                            'k2' = k_2,
                            'kd' = k_D,
                            'ke' = k_E,
                            'kg' = k_G,
                            'kr_est' = k_R_est,
                            'Total Elimination' = TotalElim_rate,
                            'kr_pct' = kr_pct,
                            'pKa' = chem$pKa,
                            'logDbw' = log10(D_BW) ,
                            'log Kow' = chem$Log_Kow,
                            'log Dmw' =log10(D_MW),
                            'log Dow' = log10(D_OW),
                            'log Kpw' = chem$Log_Kpw,
                            'D_BW' = D_BW,
                            'D_MW' = D_MW,
                            'D_OW'= D_OW,
                            'K_PW' = chem$K_PW,
                            'pHi'= env$pHi,
                            'pHg' = env$pHg,
                            'K_GB'= K_GB,
                           'nu_NB'= org$nu_NB,
                           'nu_LB'=org$nu_LB,
                           'nu_PB'=org$nu_PB,
                           'nu_OB'=org$nu_OB,
                           'nu_WB'=org$nu_WB,
                           'epsilon_N'= org$epsilon_N,
                           'epsilon_L'=org$epsilon_L,
                           'epsilon_P'=org$epsilon_P,
                           'epsilon_O'=org$epsilon_O,
                           'Log_Koc'=chem$Log_Koc,
                           'Phi'= chemdata$Phi, 
                           'RMR' = org$RMR)))
    
    return(Output_Data)
    
}

