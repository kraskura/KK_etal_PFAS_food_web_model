
# Author: Krista Kraskura
# date: feb 7 2024

# libraries 
library(tictoc)
library(here)

tic("run time for all analysis")

# 1. output the tables:
source(knitr::purl(here("Code", "R results", "model_input_data_tables.Rmd"), quiet=TRUE))

# 2. generate the food web figures and stable isotope figs
source(knitr::purl(here("Code/R results/stable_isotope_data.Rmd"), quiet=TRUE))

# 3. run the bioaccumulation main model, generate output
source(knitr::purl(here("Code/R results/main_results.Rmd"), quiet=TRUE))

# 4. run the error estimate part of the model, generate output
source(knitr::purl(here("Code/R results/sensitiv_food_web_and_metabolism.Rmd.Rmd"), quiet=TRUE))

toc()