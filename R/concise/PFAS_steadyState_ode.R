  
#' Title: Steady State Equation for calculating concentration in a single organism
#' Description: calculates species and PFAS and environment specific
#' transfer rates (uptake: k1, k_D; elimination: k2, k_E, k_R (estimated), k_G, k_M)
#' 
#' UPTAKE RATES: 
#'  k1 calculated inside the organism class
#'  kD calculated here, depends on food items (food web structure) and PFAS in each food item
#' 
#' ELIMINATION RATES: 
#'  k2 calculates here, depends on k1
#'  K_E calculated inside the organism class
#'  K_G calculated inside the organism class
#'  K_M = 0
#'
#' @param settings 
#' @param PFAA Specific PFAS analyte 
#' @param chemdata chem data Class (in 'PFAS_classes_ode.R')
#' @param env  environment Class (in 'PFAS_classes_ode.R')
#' @param chem  chem Class (in 'PFAS_classes_ode.R')
#' @param org  organism Class (in 'PFAS_classes_ode.R')
#' @param C_D  the PFAS concentrations in each food item (string of numbers, must be the same length as P)
#' @param P fraction of each food item eaten (string of numbers, must be the same length as C_D)
#' @param Pd  fraction of sediment eaten (one number)
#' @param kRTable # estimated renal elimination constants
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
SSC_B<-function(settings,
                # PFAA,
                chemdata,
                env,
                chem,
                org, 
                C_D, 
                P, 
                Pd,
                kRTable=NULL){
  
  # print("SS: new")
  
    ########################################################### *
    ## Get uptake and elimination rates -----
    ########################################################### *
    C_WDP <- chemdata$C_WDP # Freely dissolved chemical concentration in the sediment associated pore
    
    k_1 <- org$k_1 #k_1 is the uptake rate constant (L/kg*d) for chemical uptake via respiratory area
    k_2 <- org$k_2 #k_1 is the elimination rate constant (d-1) for chemical uptake via respiratory area

    if(org$switchk_1 == 0){ 
      # Phytoplankton
        k_D = 0
        k_E = 0
        k_G = org$GRF
    } else if (org$switchk_1 == 1 | org$switchk_1 == 2){ 
      # Zooplankton, Aquatic and terrestrial Invertebrates, Fish, and birds
        k_D = org$k_D
        k_E = org$k_E  #k_E is the rate constant (1/d) for chemical elimination via excretion into egested feces
        k_G = org$k_G  #k_G is the growth rate constant

    }else{
        stop('SS: error in k1 switch selection for dietary & growth uptake & elimination pathways')
    }

    if (org$switchk_R == 1) {
        # df_krkbRatio
        if(settings$chooseModel == "Liang_etal_2022"){
          k_R_est = 0
        }else{
          k_R_est = kRTable[kRTable$chemID == chem$chemID,'kr/kb'] * k_2
        }
        
    } else if (org$switchk_R == 0) {
      # plants
        k_R_est = 0
    } else if(org$switchk_R == 2){ 
      # birds
        k_R_est = 0
    }else{
      stop('SS: error in kR switch')
    }
    
    k_M <- org$k_M # metabolic transformation currently set to zero
    
    # print(c("rate constants", k_1, k_2, k_D, k_E, k_G, k_M, k_R_est))
     ############################################### *
     ## Calculate tissue concentration -----
     ############################################### *
    
    # chemdata$C_WTO in ng/mL
    # C_D ng/g
    # C_s 
    # print(c("Cw", chemdata$C_WTO))
    # print(as.numeric(unlist(C_D)))
    # print(c("Cs", chemdata$C_s))
    
    if (settings$chooseModel == "Sun_etal_2022"){
      C_B = (( (k_1 * (org$m_O * chemdata$Phi * chemdata$C_WTO + (1 - org$m_O) * C_WDP))  +
              (k_D * (sum(P * as.numeric(unlist(C_D))) + (Pd * chemdata$C_s))) ) /
              (k_2 + k_E + k_M + k_G + k_R_est))

    }else if(settings$chooseModel == "Liang_etal_2022"){
      C_B = (k_1 * chemdata$C_WTO) + 
            (k_D * (sum(P * as.numeric(unlist(C_D)))))
      
      C_B = C_B - (k_2 + k_E + k_M + k_G) * C_B
      
      message("C_B estimated using Liang et al model")
      
    }else if(settings$chooseModel == "terrestrial"){
    
      # bird no gill uptake
      C_B = (k_D * (sum(P * as.numeric(unlist(C_D))))) / (k_2 + k_E + k_M + k_G + k_R_est)
    
    }else{
      message("SS: expect error b/c model type")
    }
    
    ################################################ *
    ## Create table of output values ------
    ################################################ *
    D_OW = chem$D_OW
    D_MW = chem$D_MW
    D_BW = org$D_BW

    G_F = org$G_F # plant  = NA
    G_D = org$G_D# plant  = NaN
    G_V = org$G_V# plant  = NA
    E_D = chem$E_D
    K_GB = org$K_GB# plant  = NA
    E_W = chem$E_W

    # print(c("D_OW", as.numeric(D_OW)))
    # print(c("D_MW",as.numeric(D_MW)))
    # print(c("D_BW",as.numeric(D_BW )))
    # print(c("G_F",as.numeric(G_F )))
    # print(c("G_D",as.numeric(G_D )))
    # print(c("G_V",as.numeric(G_V )))
    # print(c("E_D",as.numeric(E_D )))
    # print(c("K_GB",as.numeric(K_GB )))
    # print(c("E_W",as.numeric(E_W )))
    
    # *************** Organism specific chemical uptake here ***********
    # if plant skip
    if(settings$chooseModel == "terrestrial"){
      Diet = sum(P * t(C_D)) # only food pfas 
      Sediment = 0 # only sediment
      Water = 0
      FeedRate = G_D / org$W_B # kg/day per kg animal
      Gill_uptake = 0
      Dietary_uptake = k_D * (sum(P * t(C_D))) # kg food/kg org * ng/g (or ug/kg food) = ug chemical/kg fish/day
      Sediment_uptake = 0 # kg food/kg org * ng/g (or ug/kg food) = ug chemical/kg fish/day
      Uptake = Gill_uptake + Dietary_uptake + Sediment_uptake
      Gill_uppct = Gill_uptake / Uptake
      Diet_uppct = Dietary_uptake / Uptake
      Sediment_uppct = Sediment_uptake / Uptake
      TotalElim_rate = (k_2 + k_E + k_M + k_G + k_R_est) 

      NL_pct = (org$nu_NB * D_OW) / D_BW
      PL_pct = (org$nu_LB * D_MW) / D_BW
      Protein_pct = (org$nu_PB * chem$K_PW) / D_BW # protein here 
      NLOM_pct = (org$nu_OB * D_OW * 0.05) / D_BW
      water_pct = (org$nu_WB) / D_BW
      
      kr_pct = k_R_est / (k_2 + k_E + k_G + k_R_est)
      
    }else{
      # if(settings$chooseMo)
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
      Dietary_uptake = k_D * (sum(P * t(C_D))) # kg food/kg org * ng/g (or ug/kg food) = ug chemical/kg fish/day
      Sediment_uptake = k_D * (Pd * chemdata$C_s) # kg food/kg org * ng/g (or ug/kg food) = ug chemical/kg fish/day
        
      Uptake = Gill_uptake + Dietary_uptake + Sediment_uptake
      Gill_uppct = Gill_uptake / Uptake
      Diet_uppct = Dietary_uptake / Uptake
      Sediment_uppct = Sediment_uptake / Uptake
      TotalElim_rate = (k_2 + k_E + k_M + k_G + k_R_est) 
      
      # *************** Organism specific chemical uptake / partitioning ***********
      NL_pct = (org$nu_NB * D_OW) / D_BW
      PL_pct = (org$nu_LB * D_MW) / D_BW
      Protein_pct = (org$nu_PB * chem$K_PW) / D_BW # protein here 
      NLOM_pct = (org$nu_OB * D_OW * 0.05) / D_BW
      water_pct = (org$nu_WB) / D_BW
  
      kr_pct = k_R_est / (k_2 + k_E + k_G + k_R_est)
    }

    # *************** Organism output table ***********
    Output_Data = t(as.data.frame(c("C_B" = C_B,
                                  'mO' = org$m_O,
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
                            'km' = k_M,
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
    # print("SS: done")
    return(Output_Data)
    
}

