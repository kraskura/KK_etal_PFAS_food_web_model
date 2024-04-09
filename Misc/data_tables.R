## li {
##   line-height: 1;
## }

## ----setup, include=FALSE-------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----source files, echo = FALSE, message=FALSE----------------------------------------------------------

library(gridExtra)
library(kableExtra)
library(here)
# here::set_here(path = "../../")

source(here("Code/R model/data_tables.R"))

source(here("Code/R data tables/WG.R"))
data.wg<-inputFiles_list
source(here("Code/R data tables/JBA_spring.R"))
data.jba.sp<-inputFiles_list
source(here("Code/R data tables/JBA_summer.R"))
data.jba.sum<-inputFiles_list



## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.wg[[2]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")


## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.jba.sp[[2]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")


## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.jba.sum[[2]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")



## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------

data.wg[[3]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")



## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------

data.jba.sp[[3]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")



## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------

data.jba.sum[[3]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")



## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.wg[[4]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")



## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.jba.sp[[4]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")



## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.jba.sum[[4]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")



## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.wg[[5]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")



## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.jba.sp[[5]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")



## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.jba.sum[[5]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")



## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.wg[[6]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")



## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.jba.sp[[6]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")



## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.jba.sum[[6]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")



## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.wg[[7]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")


## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.jba.sp[[7]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")


## ----echo = FALSE, warning=FALSE------------------------------------------------------------------------
data.jba.sum[[7]] %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Verdana")


