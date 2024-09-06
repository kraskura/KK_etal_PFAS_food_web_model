# Author: Krista Kraskura
# Date: july 17 2024
# Purpose: functions and script to transform data from eurofins format to 
# usable csv file. 

# set your own current working directory
setwd("/Users/kristakraskura/Github_repositories/KK_etal_PFAS_food_web_model/Data/SERDP_report_support/")

# libraries
library(readxl)
# library(here)
# here::i_am(path = "./Code/SERDP_report_support/data_wrangling_Eurofins_data.R")


# ******* functions ***************
# Helper function for excel names 
generate_excel_columns <- function(n_cols) {
  columns <- character(n_cols)  # Pre-allocate vector for efficiency
  
  for (i in 1:n_cols) {  # ADV is the 802nd column
    n <- i - 1  # Adjust for 0-based indexing
    column <- ""
    
    while (n >= 0) {
      column <- paste0(LETTERS[(n %% 26) + 1], column)
      n <- floor(n / 26) - 1
    }
    
    columns[i] <- column
  }
  
  return(columns)
}
  
# function to wrangle
eurofins_data_wrangle<-function(data, n_cols){
  
  contamin<-read_excel(data, range = c("A14:B54"))
  colnames(contamin)<-c("Analyte", "CAS")
  
  rows_range_start<-seq(3, n_cols, 5)
  rows_range_end<-seq(7, n_cols, 5)  
  
  for(i in 1:length(rows_range_end)){
    # the data 
    start_range<-paste(generate_excel_columns(n_cols)[rows_range_start[i]], "14", sep = "")
    end_range<-paste(generate_excel_columns(n_cols)[rows_range_end[i]], "54", sep = "")
    sample_range<-(paste(start_range, ":", end_range, sep = ""))
    sample_data<-read_excel(data, range = c(sample_range))
    # sample info, cell in excel
    sampleID_cell<-paste(generate_excel_columns(n_cols)[rows_range_start[i]], "7", sep = "")
    samplingDate_cell<-paste(generate_excel_columns(n_cols)[rows_range_start[i]], "9", sep = "")
    matrix_cell<-paste(generate_excel_columns(n_cols)[rows_range_start[i]], "10", sep = "")
    dilutionFact_cell<-paste(generate_excel_columns(n_cols)[rows_range_start[i]], "11", sep = "")
    units_cell<-paste(generate_excel_columns(n_cols)[rows_range_start[i]], "12", sep = "")
    # get excel values
    sampleID<-names(read_excel(data, range = c(sampleID_cell)))
    samplingDate<-names(read_excel(data, range = c(samplingDate_cell)))
    matrix<-names(read_excel(data, range = c(matrix_cell)))
    dilutionFact<-names(read_excel(data, range = c(dilutionFact_cell)))
    units<-names(read_excel(data, range = c(units_cell)))

    # populating sample information
    sample_data$SampleID<-sampleID
    sample_data$SamplingDate<-samplingDate
    sample_data$Matrix<-matrix
    sample_data$DilutionFact<-dilutionFact
    sample_data$Units<-units
    
    data0<-cbind(sample_data, contamin)
    data0<-data0[,c("SampleID", "SamplingDate", "Matrix", "DilutionFact", 
                    "Units", "Analyte", "CAS", "Result", "Q", "LOQ",
                    "LOD", "DL")]
    
    colnames(data0)<- c("Sample_ID", "Sampling_Date", "Matrix", "Dilution_Factor",
                        "Units", "Analyte", "CAS", "Result Value", "Q", "LOQ",
                        "LOD", "DL")
    # print(i)
    if(i == 1){
      data_new<-data0
    }else{
      data_new<-rbind(data_new,data0)
    }
  }
  
  return(data_new)
  
}
# **********************


# ******* testers ***************
generate_excel_columns(n_cols = 802)# for WG water, the last column  = ADV
generate_excel_columns(n_cols = 357)# for JBA sed, the last column  = MS
generate_excel_columns(n_cols = 837)# for JBA water, the last column  = AFE
# data<-"410-169030-1_Std_Lanc_ExcelSheet.xlsx"
# n_cols <- 802
# newdata<-read_excel("410-173731-1_Std_Lanc_ExcelSheet.xlsx")
# **********************

# **** reformatting and saving new files********
# applying the functions to convert data:
wg_water_data<-eurofins_data_wrangle(data = "410-169030-1_Std_Lanc_ExcelSheet.xlsx", 
                                     n_cols = 802) # 160 samples (5*160)+2 = 802

sed_data<-eurofins_data_wrangle(data = "410-173731-1_Std_Lanc_ExcelSheet.xlsx", 
                                     n_cols = 357) # 71 samples (5*71)+2 = 357

jba_water_data<-eurofins_data_wrangle(data = "410-165666-1_Std_Lanc_ExcelSheet.xlsx", 
                                     n_cols = 837) # 167 samples (5*167)+2 = 837

jba_sed_data<-sed_data[grepl("JBA", sed_data$Sample_ID), ]
wg_sed_data<-sed_data[grepl("WG", sed_data$Sample_ID), ]
  
# save files 
write.csv(wg_water_data, file = "410-169030_WG_water_data_reformatted.csv", row.names = FALSE)
write.csv(jba_water_data, file = "410-165666_JBA_water_data_reformatted.csv", row.names = FALSE)
write.csv(wg_sed_data, file = "410-173731_WG_sediment_data_reformatted.csv", row.names = FALSE)
write.csv(jba_sed_data, file = "410-173731_JBA_sediment_data_reformatted.csv", row.names = FALSE)
