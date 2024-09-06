# **PFAS Bioaccumulation in fish.**

Materials used for SERDP report: September 6, 2024.

Materials published: Part of the data are published in [Brown et al 2023](https://doi.org/10.1016/j.scitotenv.2023.163149)

Authors: Krista Kraskura (kkraskura\@towson.edu), Abbi Brown, Christopher Salice.

------------------------------------------------------------------------

### **PFAS food web bioaccumulation model**

The adopted PFAS food web was published by [Sun et al. (2022)](https://doi.org/10.1039/D2EM00047D) and was forked from <https://github.com/SunderlandLab/fish_foodweb_pfas_model>. The model was reprogrammed in R and implemented using data collected from various field studies.

The model use and implementation: The code **EXECUTE.R** in the local directly runs all modelling files to obtain results.

##### .Code/R model/

-   contains R files (all functions):

    -   "PFAS_steadyState_v2.R": a script/function that estimates the steady state chemical specific PFAS uptake. The uptake from sediment, water, and diet are estimated using this function. Focuses on PFAS chemical at the time.

    -   "PFAS_bioaccum_v2.R": a script/function with bioaccumulation model for each species (integrates the food web structure). The steady state model is applied to each species. Focuses on PFAS chemical at the time.

    -   "runEcosystemModels.R": a script/function that runs the bioaccumulation model across the following scenarios, median or mean PFAS conditions, min and max PFAS value conditions. Runs a loop to cover all PFAS as provided by the input data.

    -   "PFAS_classes_v2.R": a script (set of functions) where all model classes get defined. These set of functions define environment, chemical, and organism. Several key equations outlining organismal biology are defined here (e.g., feeding rate, ventilation rate). The rest of the model depends on these classes.

    -   "data_tables.R": a script/function that formats the data tables (input files in the model)

    -   "data_MonteCarlo.R": a script/function that consolidates and prepares data for Monte Carlo simulations.

    -   "runErrorEstimates.R": a script/function that estimates average model bias, coefficient of model bias, the R2 at the whole model level, or PFAS and species specific.

    -   "mc_sim_data_table.R": a script/function that runs a mc simulation with input data from data_MonteCarlo.R.

##### .Code/R Results/

-   R markdown files where all model results are obtained. The scripts rely on the defined functions in Code/R model.

##### .Code/R data tables/

-   Script to construct all input data tables.

    -   Values are first extracted from documented data files (found in /Data/Original_Brown_etal/ folder

    -   Then the measured environmental values, along with chemistry specific values (e.g., Koc), and food web are provided in the function that is defined earlier.

    -   These input data tables are used to parameterize the model.

##### .Data/

-   The orginal data files with all measured fish, sediment, and water PFAS values.

-   Additional, supporting data files to execute the code (e.g., the renal elimination factors reported and published by Sun et al 2022).

-   Output files from the Bioaccumulation part (part 2, see below), the outputs for the model were too large and not included with the GitHub version control.

##### .Figure/

-   All saved output figures.

------------------------------------------------------------------------

### **Fish bioaccumulation analysis**

-   Input data (original data) as well as few output files (statistic summaries) for the bioaccumulation part are found in .Data/SERDP_report_support/.

-   R scripts to support creation of reported figures and provide all analyses are in .Code/SERDP_report_support

-   The output figures are found in ./Figure/SERDP/

------------------------------------------------------------------------

### Contact info:

Please contact Krista Kraskura at kkraskura\@towson.edu with any questions and concerns.

Please contact us if interested in using any data or script published in this repository.
