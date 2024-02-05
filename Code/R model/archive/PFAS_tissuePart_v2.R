# EXTRA CALCULATIONS -----
# Extract partitioning coefficients
choose_partitioning <- function(PFAS_list,  inputFiles_list){

    chemicalParams <- as.data.frame(inputFiles_list$chemicalParams[PFAS_list])
    Log_Kow <- chemicalParams['Log_Kow',]
    pKa <- chemicalParams['pKa',]
    Log_Dmw <- chemicalParams['log_Dmw',]
    Log_Kpw <- log10(1.36 * 10^chemicalParams['log_Kpw',])
    K_OW <- 10^Log_Kow #added from the class definitions; getD_OW funct

    # calculate Dow
    getD_OW <- function(pHi, pKa, K_OW){
      Xn = 1 / (1 + 10^(pHi-pKa)) # neutral fraction of the compound
      Xi = 1 - Xn # ionic fraction of the compound

      logK_OW_i <- log10(K_OW) - 3.1 # based on Armitage et al 2013, see SI Table S7
      K_OW_i <- 10^logK_OW_i
      D_OW = Xn * K_OW + Xi * K_OW_i

      if(log10(D_OW) < 0){
        message("error in D_OW calculation")
      } else {
        return(D_OW)
      }
    }

    pHi <- 7.4


    for (i in 1:length(PFAS_list)){
        Dow <- getD_OW(pHi, pKa[PFAS_list[i]], 10^Log_Kow[PFAS_list[i]])
        if(i == 1){
          Log_Dow <- c(log10(Dow))
        } else {
          Log_Dow <- append(Log_Dow, log10(Dow))
        }
    }

    # Log_Dow <- pd.Series(Log_Dow, index=PFAS_list)
    partitionCoeff <- data.frame('logKow' = t(Log_Kow), 'logDow' = t(as.data.frame(Log_Dow)),
                                 'logDmw' = t(Log_Dmw), 'logKpw' = t(Log_Kpw))
    names(partitionCoeff) <- c('logKow', 'logDow',
                                'logDmw', 'logKpw')
    return(t(partitionCoeff))
}


# Calculate tissue contributions
calc_TissueContributions <- function(PFAS_list,  inputFiles_list, Spp){
    tissuePct <- inputFiles_list$organismData

    DataTable <- choose_partitioning(PFAS_list, inputFiles_list = inputFiles_list)

    Total_Dbw<- (10^DataTable['logDow',] * tissuePct['nu_NB', Spp]) +
      (10^DataTable['logDmw',] * tissuePct['nu_LB', Spp]) +
      (10^DataTable['logKpw',] * tissuePct['nu_PB', Spp]) +
      (10^DataTable['logDow',] * 0.05 * tissuePct['nu_OB',Spp]) +
      tissuePct['nu_WB',Spp]

    Total_log_Dbw <- log10(Total_Dbw)

    Lipid_p <- (10^DataTable['logDow',] * tissuePct['nu_NB',Spp]) / Total_Dbw
    Phospholipid_p <- (10^DataTable['logDmw',] * tissuePct['nu_LB',Spp]) / Total_Dbw
    Albumin_p <- (10^DataTable['logKpw',] * tissuePct['nu_PB',Spp]) / Total_Dbw
    NLOM_p <- (10^DataTable['logDow',] * 0.05 * tissuePct['nu_OB',Spp]) / Total_Dbw
    Water_p <- tissuePct['nu_WB',Spp] / Total_Dbw

    DataTable<-rbind('Total (Dbw)' = t(data.frame(Total_Dbw)),
                        'Total (log Dbw)' = t(data.frame(Total_log_Dbw)),
                        'Lipid%' = t(data.frame(Lipid_p)),
                        'Phosphoipid%' = t(data.frame(Phospholipid_p)),
                        'Albumin%' = t(data.frame(Albumin_p)),
                        'NLOM' = t(data.frame(NLOM_p)),
                        'Water%' = t(data.frame(Water_p)),
                    DataTable)
    rownames(DataTable)<-c("Total (Dbw)", "Total (log Dbw)",
                           "Lipid%","Phospholipid%","Albumin%",
                           "NLOM%", "Water%",
                           "logKow", "logDow", "logDmw" ,"logKpw")

    return(invisible(DataTable))
}
