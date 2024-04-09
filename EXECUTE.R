
# Author: Krista Kraskura
# date: feb 7 2024

# libraries 
library(tictoc)
library(here)

tic("run time for all analysis")

# 1. output the tables:
source(knitr::purl(here("Code", "R results", "data_tables.Rmd"), quiet=TRUE))

# 2. generate the food web figures and stable isotope figs
source(knitr::purl(here("Code/R results/food_web.Rmd"), quiet=TRUE))

# 3. run the bioaccumulation main model, generate output
source(knitr::purl(here("Code/R results/objectives_results.Rmd"), quiet=TRUE))

# 4. run the error estimate part of the model, generate output
source(knitr::purl(here("Code/R results/objectives_results3_errorsEst.Rmd"), quiet=TRUE))

toc()