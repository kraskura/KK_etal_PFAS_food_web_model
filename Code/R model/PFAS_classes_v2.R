
# Author: Krista Kraskura
# Date: march 2024
# Description: Creating new data classes that define the physical and biological
# environments, the organism, the chemical (different PFAS)
# Use:


# # install.packages('assertthat')
library(assertthat)

# Template to create classes:

# new_Class<-function(x, env){
#  	# validation 
#   assert_that(all(!is.na(as.numeric())))# 	
#   # attributes 
# 	x = x
# 	
# 	# structure
# 	structure(
# 		list(
# 		  get_vars = function(){
# 		    list("x" = x)
# 		  }
# 	  ), 
# 		class = "new_Class"	
# 	)	
# }

 
# Data classes for input parameters and models ----
# object oriented implementations

# ***********************************************
new_Settings <- function(
		chooseDiet, # default (trophic ) or forced (known diet); 
		chooseDmw = "Droge", # or 'calculated'
		chooseEw = "empirical",
		chooseKoc = "Koc",
		chooseModel = "Sun_etal_2022"){
	
  # validation
	stopifnot(is.character(chooseDiet),
		is.character(chooseDmw),
		is.character(chooseEw),
		is.character(chooseKoc),
		is.character(chooseModel))
		
	# structure 
	structure(
		list(
		chooseDiet = chooseDiet,
		chooseDmw = chooseDmw,
		chooseEw = chooseEw,
		chooseKoc = chooseKoc,
		chooseModel = chooseModel),
		class = "new_Settings"
	)
}

# ***********************************************
new_Environment <- function(
    C_OX,
    T,
    C_SS,
    pH,
    pHg,
    OCS){
	# validation 
  assert_that(all(!is.na(as.numeric(C_OX, T, C_SS, pH, OCS))))
	
  # attributes 
	C_OX = C_OX
	T = T
	C_SS = C_SS
  OCS = OCS
	pH = pH
	pHi = 7.4 # =constant 
	pHg = pH-1 
	
	# structure
	structure(
		list(
		  get_vars = function(){
		    list("C_OX" = C_OX, "T" = T, "OCS" = OCS, "C_SS" = C_SS,
		         "pH" = pH, "pHi" = pHi, "pHg" = pHg)
		  }
	  ), 
		class = "new_Environment"	
	)	
}

# ***********************************************
new_Chemical <- function(
    Log_Kow,
    Log_Kpw,
    Log_Dmw,
    pKa,
    Ed_exp,
    Ew_exp,
    Log_Koc,
    chemID,
    env,
    settings){
  
  #   """ Dataclass that holds parameters unique to each PFAS chemical.
  #   Parameters include physiochemial properties, pH-dependent properties (varies with environment), 
  #   and empirically-derived parameter values (sometimes varies with organism characteristics/study;
  #   this choice is made during the instantiation of this class)
      
  #   Functions include calculation of chemical-specific parameters,
  #   which may depend on environmental parameters (class Environment)"""

  # validator
  assert_that(all(!is.na(as.numeric(
    Log_Kow, Log_Kpw, Log_Dmw,
    pKa, Ed_exp, Ew_exp, Log_Koc,
    chemID # in data csv provided as integers (e.g. 0, 2, 3 ..), numeric in R TRUE
  ))))
	
  # attributes 
  Log_Kow = Log_Kow 
  Log_Kpw = Log_Kpw 
  Log_Dmw = Log_Dmw 
  pKa = pKa
  Ed_exp = Ed_exp
  Ew_exp = Ew_exp
  Log_Koc = Log_Koc
  chemID = chemID 
  env = env # call the environment with defined attributes
	
  K_OW = 10^Log_Kow
  K_PW = 10^Log_Kpw
  D_MW_exp = 10^Log_Dmw
  K_OC = 10^Log_Koc
    
  # def getD_OW(self, env: Environment) -----
  #   """ Octanol-water distribution coefficient (Armitage et al. 2013)""" 
  Xn = 1 / (1 + 10^(env$pHi-pKa)) # neutral fraction of the compound
  Xi = 1 - Xn # ionic fraction of the compound
  logK_OW_i = log10(K_OW) - 3.1 # based on Armitage et al 2013, see SI Table S7 # =constant
  K_OW_i = 10^logK_OW_i
  D_OW = Xn * K_OW + Xi * K_OW_i # =return

  assert_that(log10(D_OW) > 0, msg = 'error in D_OW calculation (negative value)')
  
  # def getD_MW(self, env: Environment, settings: Settings): --------
  if (settings$chooseDmw == 'Droge'){
      D_MW = D_MW_exp          
  } else if (settings$chooseDmw == 'calculated'){
    
      Xn = 1 / (1 + 10^(env$pHi-pKa)) # neutral fraction of the compound
      Xi = 1 - Xn # ionic fraction of the compound

      a = 1.01 # Endo et al, as cited in Armitage et al 2013, Table S2 # =constant
      b = 0.12 # Endo et al, as cited in Armitage et al 2013 # =constant
      logK_MW_n = (a * log10(K_OW)) + b
      K_MW_n = 10^logK_MW_n

      logK_MW_i = log10(K_MW_n) - 2 # based on Armitage et al 2013, see SI Table S7
      K_MW_i = 10^logK_MW_i

      D_MW = Xn * K_MW_n + Xi * K_MW_i # =return
      
  } else {
      stop('error in D_MW choice')
  }
  assert_that(log10(D_MW) > 0, msg = "error in D_MW calculation (negative value)")

    
  # def getE_W(self, env: Environment, settings: Settings): --------
  #     """ Fish gill uptake efficiency (%)
  #     Used to calculate k1 clearance rate constant for fish, invertebrates, and zooplankton
  #     """
  if(grepl('beta', settings$chooseEw)){
      pHg = env$pHi - 1
      Ew_mu = 1 + 10^(pHg - pKa)
      Ew_beta = as.numeric(gsub('beta','', settings$chooseEw)) 
      assert_that(is.numeric(Ew_beta), msg = "Ew_beta not numeric.
                  new_Settings$chooseEw must be in format: beta#, #beta
                  [# = numeric var, no spaces, no special char]")
      if(settings$chooseModel == "Sun_etal_2022"){
        E_W = (1 / (1.85 + (Ew_mu * 155/K_OW))) + Ew_beta # =return
        
      }else if(settings$chooseModel == "Liang_etal_2022"){
        E_W = (1 / (1.85 + (155/K_OW))) 
      }else{
        message("expect error b/c chooseModel choice")
      }
      
  } else if (settings$chooseEw == 'empirical'){
      E_W = Ew_exp # =return
  } else {
      stop('error in Ew selection')
  }

  # def getE_D(self, settings: Settings):  ----------
  #     """ Dietary chemical transfer efficiency (%)
  #     Used to calculate kd uptake rate constant.
  #     Empirical values for PFAAs differ for juvenile and adult fish
  #     """
  if(settings$chooseModel == "Sun_etal_2022"){
    E_D = Ed_exp # =return
  }else if(settings$chooseModel == "Liang_etal_2022"){
    E_D = 1 / ((5.1 * 10^(-8) * K_OW + 2))
  }else{
    message("expect error b/c model type")
  }

	structure(
		list(
		  get_vars = function(){
		    list(
		      "Log_Kow" = Log_Kow, "Log_Kpw" = Log_Kpw, "Log_Dmw" = Log_Dmw,
		      "pKa" = pKa, "D_MW_exp" = D_MW_exp, "Ed_exp" = Ed_exp,
		      "Ew_exp" = Ew_exp, "Log_Koc" = Log_Koc, "chemID" = chemID,
		      "D_OW" = D_OW, "D_MW" = D_MW, "E_W" = E_W, "E_D" = E_D,
          "K_OW" = K_OW, "K_PW" = K_PW, "K_OC" = K_OC
		      )
		  }
	  ), 
		class = "new_Chemical"	
	)	
}


# ***********************************************
new_Organism <- function(
    W_B, m_O, GRF, nu_NB, nu_LB, nu_PB, nu_OB,
    nu_WB, nu_ND, nu_LD, nu_PD, nu_OD, nu_WD,
    epsilon_N, epsilon_L, epsilon_P, epsilon_O, epsilon_W,
    switchk_1,
    switchG_D,
    switchk_R,
    A, B,
    RMR,
    MMR, 
    P_B,
    feeding_exp,
    scaling_exp,
    sigma,
    env, settings, chem){
  # """ Dataclass that holds parameters unique to each organism, or the calculation of parameters that may vary
  # depending on organism type (i.e. sigma, beta)
  # 
  # Functions include calculations of 
  # (1) diet composition; 
  # (2) bioenergetic rates;
  # (3) partitioning coefficients;
  # (4) chemical uptake and elimination rate constants"""
  
  # validation 
  assert_that(all(!is.na(as.numeric(
    m_O, GRF, nu_NB, nu_LB, nu_PB, nu_OB,
    nu_WB, nu_ND, nu_LD, nu_PD, nu_OD, nu_WD,
    epsilon_N, epsilon_L, epsilon_P, epsilon_O, epsilon_W,
    switchk_1, switchG_D, switchk_R, # integers 
    A, B, feeding_exp, scaling_exp, sigma))),
    msg = "new_Organism: input vars not numeric")
  assert_that(all(c(is.na(W_B) | !is.na(is.numeric(W_B)))),
    msg = "new_Organism: input vars not numeric") # phytoplankton can be body weight NA
	
  # attributes 
  W_B = W_B
  m_O = m_O
  GRF = GRF
  
  scaling_exp = scaling_exp
  feeding_exp = feeding_exp
  MMR = MMR
  RMR = RMR

  nu_NB = nu_NB
  nu_LB = nu_LB
  nu_PB = nu_PB
  nu_OB = nu_OB
  nu_WB = nu_WB

  nu_ND = nu_ND
  nu_LD = nu_LD
  nu_PD = nu_PD
  nu_OD = nu_OD
  nu_WD = nu_WD

  epsilon_N = epsilon_N
  epsilon_L = epsilon_L
  epsilon_P = epsilon_P
  epsilon_O = epsilon_O
  epsilon_W = epsilon_W

  switchk_1 = switchk_1
  switchG_D = switchG_D
  switchk_R = switchk_R

  A = A
  B = B
  sigma = sigma
  P_B = P_B
	
  # Diet functions------ 
  # def _getnu_NG(self):------ 
  nu_NG = (1 - epsilon_N) * nu_ND / ((1- epsilon_N) * nu_ND + (1 - epsilon_L) * nu_LD +
                              (1 - epsilon_P)* nu_PD + (1 - epsilon_O) * nu_OD + (1 - epsilon_W) * nu_WD)

  # def _getnu_LG(self): ------ 
  nu_LG = (1 - epsilon_L) * nu_LD / ((1- epsilon_N) * nu_ND + (1 - epsilon_L) * nu_LD +
                                  (1 - epsilon_P)* nu_PD + (1 - epsilon_O) * nu_OD + (1 - epsilon_W) * nu_WD)

  # def _getnu_PG(self): ------
  nu_PG = (1 - epsilon_P) * nu_PD / ((1- epsilon_N) * nu_ND + (1 - epsilon_L) * nu_LD +
                                  (1 - epsilon_P)* nu_PD + (1 - epsilon_O) * nu_OD + (1 - epsilon_W) * nu_WD)

  # def _getnu_OG(self): ------ 
  nu_OG = (1 - epsilon_O) * nu_OD / ((1- epsilon_N) * nu_ND + (1 - epsilon_L) * nu_LD +
                                  (1 - epsilon_P)* nu_PD + (1 - epsilon_O) * nu_OD + (1 - epsilon_W) * nu_WD)

  # def _getnu_WG(self): ------ 
  nu_WG = (1 - epsilon_W) * nu_WD / ((1- epsilon_N) * nu_ND + (1 - epsilon_L) * nu_LD +
                                  (1 - epsilon_P)* nu_PD + (1 - epsilon_O) * nu_OD + (1 - epsilon_W) * nu_WD)
      
  # = return: nu_NG, nu_LG nu_OG, nu_NG, nu_WG

  # Organism functions --------
  # def _getG_V(self, env: Environment, settings: Settings): ----
  #     "G_V is the ventilation rate. C_OX is the dissolved oxygen concentration (mg/L)"
  
  # G_V: ventilation rate = 1400 * (W_B^0.65) / env$C_OX" (unit: L/d)
  # Vox: metabolic rate = 980 * WB^0.65 ;units: mgO2/day from Arnot and Gobas 2004
  if(!is.na(RMR)){    # RMR in units mg O2 / day 
    if (RMR == 0){

      T = env$T
      
      # Robinson et al (1983) The effects of body size and temperature on metabolic rate of organisms
      RMR_ml0 = (0.067 * (W_B*1000)^(scaling_exp-1) * exp(0.051*T)) # mlO2 /g /h
      RMR_ml = RMR_ml0 * (W_B*1000) * 24 # mlO2 /whole fish / day
      RMR_n = ((1) * (RMR_ml / 1000)) / ((0.0821) * (273+T)) # mol O2
      RMR =  RMR_n * 32* 1000 # mg O2 / whole fish / day
      
      # P*V = n*R*T / n = 
      # where P = pressure in atm (1 atm = 1000mb = 101.3 kPa), 
      # V = vol in Liters
      # n = moles of gas
      # gas consatnt (for these units, =0.0821 L*atm/mol*K)
      # temperature in Kelvins ( T(K) = T + 273K )
      # mol weight O2 32.0 g/mol
      # 
      # RMR0 = (W_B*1000)^(scaling_exp-1) * exp(0.06*T)  # g O2 / g/ d # fish bioenergetics 4.0
      # RMR = RMR0 * 1000 
      
      G_V = RMR / (env$C_OX * m_O)
      
    }else{
      # message("RMR provided")
      G_V = RMR/ env$C_OX * m_O
    }
    
  }
  if (is.na(RMR)){
    # message("G_V original model, no RMR")
    if(settings$chooseModel == "Sun_etal_2022"){
      G_V = 1400 * (W_B^scaling_exp) / env$C_OX  
    }else if(settings$chooseModel == "Liang_etal_2022"){
      G_V = 980 * (W_B^0.65) / env$C_OX # changed
    }else{
      message("expect error b/c model type")
    }
    

  }


  # 2. Mod 1 (see link: http://rem-main.rem.sfu.ca/papers/gobas/A%20Review%20of%20Bioconcentration%20factor%20(BCF)%20and.pdf
  # G_V = 980 * (W_B^0.65) / env$C_OX 
  
  # 3. Mod 2. (see link: https://deepblue.lib.umich.edu/bitstream/handle/2027.42/23722/0000694.pdf?sequence=1
  # G_V = 1400 * (W_B^0.8) / env$C_OX 
  
  # def _getG_D(self, env: Environment, settings: Settings): ----
  #     """ G_D is the feeding rate
  #     G_Dswitch determines which formula to use based on the type of organism and data availability
  #     1 should be used to estimate G_D for coldwater fish species, and in some cases, zooplankton
  #     2 should be used to estimate G_D for filter-feeding species
  #     Use empirical data where available,
  #     though be careful about the energetic content of food and associated growth rates"""

  if(switchG_D == 0){
    G_D = NaN
  } else if (switchG_D == 1){
  # Note that several studies have observed decreasing concentrations with increasing fish length w/in the same species
  # This could potentially be due to growth dilution, onotogenetic shifts, or other processes, not necessarily feeding rate

    if(settings$chooseModel == "Sun_etal_2022"){
      G_D = 0.022 * W_B^feeding_exp * exp(0.08*env$T) # McLeod et al 2016 for 0.08*T 
    }else if(settings$chooseModel == "Liang_etal_2022"){
      G_D = 0.22 * W_B^0.85 * exp(0.06*env$T) # changed
    }else{
      message("expect errro b/c model type")
    }

    # return G_D
  } else if (switchG_D == 2){
      #G_V is the gill ventilation rate
      #C_SS is the concentration of suspended solids
      #sigma is the scavenging efficiency of particles
      # G_V = _getG_V(env, settings)
      G_D = G_V * env$C_SS * sigma # return G_D
      message("G_D switch = 2")
  } else {
      stop("error in getG_D: wrong input for switchG_D argument")
  }
  
  # def _getG_F(self, env: Environment, settings: Settings):  -----
  # G_F is the fecal egestion rate
  # epsilon_L, epsilon_N, and epsilon_W are the dietary assimilation efficiencies 
  # of lipid, NLOM and water, respectively.
  # nu_LD, nu_ND, and nu_WD are the overall lipid,
  # NLOM and water contents of the diet, respectively.
  # G_D = self._getG_D(env, settings) # defined above
  G_F = ((1 - epsilon_N) * nu_ND +
           (1 - epsilon_L) * nu_LD +
           (1 - epsilon_P) * nu_PD +
           (1 - epsilon_O) * nu_OD +
           (1 - epsilon_W) * nu_WD) * G_D # return G_F

  # Partitioning coefficients ---------
  # def getD_BW (self, chem: Chemical, env: Environment, settings: Settings):------ 
  D_OW = chem$D_OW
  D_MW = chem$D_MW

  if(switchk_1 == 0){
    Beta = 0.35
  } else {
    Beta = 0.05
  }
  
  D_BW = nu_NB * D_OW + nu_LB * D_MW + nu_PB * chem$K_PW + nu_OB * D_OW * Beta + nu_WB # = return
 
  # def _getK_GB(self, chem: Chemical, env: Environment, settings: Settings): -----
  #     """ K_GB is the partition coefficient of the chemical between the contents of GIT and the oragnism.
  #     nu_NG, nu_LG, nu_PG and nu_WG are the neutral lipid, phospholipid, protein and water contents,respectively, in the gut."""

  
  D_GW = nu_NG * D_OW + nu_LG * D_MW + nu_PG * chem$K_PW + nu_OG * D_OW * 0.05 + nu_WG
  K_GB = D_GW / D_BW # = return

  # def _getK_BU (self, chem: Chemical, env: Environment, settings: Settings): ----
  #K_BU is the partition coefficient between the organism and urine C_B/C_U
  #K_BW is the organism water partition coefficient
  #nu_LB is the lipid fraction (kg lipid/kg organisms ww)
  #nu_NB is the NLOM fraction (kg NLOM/kg organism ww)
  #nu_WB is the water content (kg water/kg organism ww) of the organism
  #Beta is a proportionality constant expressing the sorption capacity of NLOM to that of octanol
  #.035 is a reasonable estimate of Beta

  LipidDensity=0.9 # kg/L
  PLipidDensity=0.9 # kg/L (Armitage et al. 2013)
  ProteinDensity= 1.35 # g/cm3 = kg/L (Fischer et al. 2009; Allendorf et al. 2020)
  NLOMDensity=1 # kg/L

  # below are already defined above 
  # D_OW = chem.getD_OW(env)
  # D_MW = chem.getD_MW (env, settings)

  K_BU = (nu_NB * D_OW / LipidDensity) +
    (nu_LB * D_MW / PLipidDensity) +
    (nu_PB * chem$K_PW / ProteinDensity) +
    (nu_OB * D_OW * 0.05 / NLOMDensity) + nu_WB  # return K_BU

  # Rate constant functions ---------------
# 
# def getk_1(self, chem: Chemical, env: Environment, settings: Settings): -----
#     """ k_1 is the clearance rate constant (L/kg*d) 
#         for chemical uptake via respiratory area (i.e., gills and skin) """
# 
    if(switchk_1 == 0){
      # #A and B are constants describing the resistance to chemical uptake
      # #through aqueous (A) and organic phases (B) of the algae, phytoplankton, or macrophyte.
      # #A default value = 6x10^-5 for PCBs in Arnot & Gobas 2004; 8.2 x 10^-3 in Gobas & MacKay 1987
      # #B default value = 5.5 for PCBs in Arnot & Gobas 2004; 0.68 in Gobas & MacKay 1987
      # # difference between these two options is minimal
      # 
      # #k_1 for algae, phytoplankton, and aquatic macrophytes
      k_1 = 1/(A + (B / chem$D_MW_exp)) # = return k_1

    } else if (switchk_1 == 1){
      #G_V is ventilation rate (should be able to find and not use equation)
      #G_V can be approximated by G_V = 1400 * (W_B^.65)/C_OX where C_OX = (-.24*T + 14.04) * S
      #W_B is the wet weight of the organism (kg)
      #E_W is the gill chemical uptake efficiency
      
      E_W = chem$E_W
      
      if(settings$chooseModel == "Sun_etal_2022"){
        # G_V = self._getG_V(env, settings) # available above
        #k_1 for fish, invertebrates, and zooplankton
        k_1 = E_W * G_V / W_B # =return
      }else if(settings$chooseModel == "Liang_etal_2022"){
        k_1 = (E_W * G_V) / W_B # same for both 
      }else{
        message("anticipate error, need `settings$chooseModel`")
      }
    }

# def getk_D(self, chem: Chemical, env: Environment, settings: Settings): -----
#     """ Dietary uptake clearance constant (kg food/kg organism per day)
#     The rate at which chemicals are absorbed from the diet via the GIT"""
#     #E_D is the dietary chemical transfer efficiency
#     #G_D is the feeding rate (kg food/day) - empirical data often available
#     #W_B is the weight of the organism (kg organism)
#     G_D = self._getG_D(env, settings) # available above
  
  E_D = chem$E_D
  k_D = E_D * G_D / W_B # =return k_D dietary uptake: kg/kg/d
  # both model types  the same 


# def getk_E(self, chem: Chemical, env: Environment, settings: Settings): -----
#     """ Fecal elimination rate constant: the rate at which chemicals are eliminated by the egestion of fecal matter (1/d) """
#     #G_F is the fecal egestion rate (kg-feces/kg-organism/day)
#     #E_D is the dietary chemical transfer efficiency
#     #K_GB is the partition coefficient of the chemical between the GIT and the organism.
#     #W_B is the weight of the oyrganism
#     #epsilon_L, epsilon_N, and epsilon_W are the dietary assimilation efficiencies of lipid, NLOM, and water, respectively
#     #nu_LD, nu_ND, and nu_WD are the overall lipid, NLOM, and water contents of the diet
#     G_F = self._getG_F(env, settings)
#     K_GB = self._getK_GB(chem, env, settings) 

  E_D = chem$E_D
  
  # k_e exretion 
  if(settings$chooseModel == "Sun_etal_2022"){
    k_E = G_F * E_D * K_GB / W_B # =return k_E 1/day 
  }else if(settings$chooseModel == "Liang_etal_2022"){
    k_E = 0.25 * k_D # 1/day
  }else{
    message("expect error b/c model type")
  }
  
  # k_G - feeding rate
  if (switchk_1 == 0){ # algae, phytoplankton, macrophytes
      k_G = GRF # GRF input for phytoplankton is simply the growth rate constant
  } else if (switchk_1 == 1){ # fish, invertebrates
    if(settings$chooseModel == "Sun_etal_2022"){
      k_G = GRF * W_B^-0.2 # return k_G
    }else if(settings$chooseModel == "Liang_etal_2022"){
      k_G = GRF * W_B^-0.2
    }else{
      message("error b/c model type")
    }
      
  } else {
      stop('error in getk_G: wrong input for switchk_1 argument')
  }
  
	# structure
	structure(
		list(
		  get_vars = function(){
		    list("W_B" = W_B, "m_O" = m_O, "GRF" = GRF,
		         "nu_NB" = nu_NB, "nu_LB" = nu_LB, "nu_PB" = nu_PB, "nu_OB" = nu_OB,
             "nu_WB" = nu_WB, "nu_ND" = nu_ND, "nu_LD" = nu_LD, "nu_PD" = nu_PD,
		         "nu_OD" = nu_OD, "nu_WD" = nu_WD,
		         "epsilon_N" = epsilon_N, "epsilon_L" = epsilon_L, "epsilon_P" = epsilon_P,
		         "epsilon_O" = epsilon_O, "epsilon_W" = epsilon_W,
		         "switchk_1" = switchk_1, "switchG_D" = switchG_D, "switchk_R" = switchk_R, # integers 
		         "A" = A, "B"= B, "sigma" = sigma, 
		         "nu_NG" = nu_NG, "nu_LG" = nu_LG, "nu_OG" = nu_OG, "nu_NG" = nu_NG, "nu_WG" = nu_WG,
		         "G_F" = G_F, "D_BW" = D_BW, "K_GB" = K_GB, "K_BU" = K_BU,
		         "k_1" = k_1, "k_D" = k_D, "k_E" = k_E, "k_G" = k_G, "G_D" = G_D, "G_V" = G_V, 
		         "MMR" = MMR, "RMR" = RMR, "P_B" = P_B,
		         "scaling_exp" = scaling_exp, "feeding_exp" = feeding_exp)
		  }
	  ), 
		class = "new_Organism"	
	)	
	
	
}

# test class outputs
# organism<-new_Organism(W_B = 0.5, m_O = 1, GRF = 0,
#              nu_NB = 0.04, nu_PB = 0.003, nu_LB = 0.1, 
#              nu_OB = 0.15, nu_WB = 0.8,
#              nu_ND = 0.1, nu_PD = 0.02, nu_LD = 0.3, 
#              nu_OD = 0.5, nu_WD = 2,             
#              epsilon_N = 0.9, epsilon_L = 0.9, epsilon_P = 0.9,
#              epsilon_O = 0.6, epsilon_W = 0.7,
#              switchk_1 = 1, switchG_D = 1, switchk_R = 1, A = NA, B = NA, sigma = 0.5,
#              env = env1, settings = sett, chem = chem)$get_vars()

# ***********************************************
new_ChemData <- function(
    C_WTO,
    C_s,
    Phi,
    env,
    chem){
  # """ Dataclass that holds chemical data that is unique to each model system. This includes both model parameters and input variables"""
  # #Cs is the Concentration in sediment
  # #OCS is the organic carbon content of the sediment
  # #Koc is the organic carbon partitioning coefficient; Koc = Kd / OCS 
 	
  # validation
  assert_that(all(!is.na(as.numeric(C_WTO, C_s, Phi, env, chem))))
  # attributes
	C_WTO = C_WTO
	C_s = C_s
	Phi = Phi # should be < 1; default is 1 for the lab studies if this is not used
  
	#def getC_WDP(self, chem: Chemical, env: Environment, settings: Settings): ----
  # """ Freely dissolved chemical concentration in the sediment associated pore (or interstitial) water (g/L)"""
  # C_WDP = C_SOC / K_OC where C_SOC is the chemical concentration in sediment normalized for OC content (g/kg-OC)

  C_WDP = (C_s / env$OCS) / chem$K_OC

	# structure
	structure(
		list(
		  get_vars = function(){
		    list("C_WTO" = C_WTO, "C_s" = C_s, "Phi" = Phi, "C_WDP" = C_WDP)
		  }
	  ),
		class = "new_ChemData"
	)
}

# test new_ChemData outputs--------
# chemDat<-new_ChemData(C_WTO = 0.0005, C_s = 0.5, Phi = 0.5, env = env1, chem = chem)$get_vars()



