## ----setup, cache = FALSE, echo = FALSE, message = FALSE, warning = FALSE, tidy = FALSE----
knitr::opts_chunk$set(eval = FALSE)

## ----gravity_no_zeros----------------------------------------------------
#  # 1: Import and read the dataset
#  
#  url <- "http://econ.sciences-po.fr/sites/default/files/file/tmayer/data/col_regfile09.zip"
#  zip <- "col_regfile09.zip"
#  
#  if (!file.exists(zip)) { try(download.file(url, zip)) }
#  try(system("7z e -aos col_regfile09.zip"))
#  
#  library(haven)
#  col_regfile09 <- read_dta("col_regfile09.dta")
#  
#  # 2: Isolation of one year
#  
#  library(dplyr)
#  data06 <- col_regfile09 %>%
#    filter(year == 2006)
#  
#  # 3: Choosing variables
#  
#  data06 <- data06 %>%
#    select(iso_o, iso_d, distw, gdp_o, gdp_d, rta, flow, contig, comlang_off, comcur)
#  
#  # 4: Isolation of complete cases
#  
#  library(tidyr)
#  gravity_zeros <- data06 %>%
#    drop_na()
#  
#  # 5: Exclusion of trade flows equal to 0
#  
#  gravity_no_zeros <- gravity_zeros %>%
#    filter(flow > 0)
#  
#  # 6: Export the data
#  
#  save(gravity_zeros, file = "gravity_zeros.rdata", compress = "xz")
#  save(gravity_no_zeros, file = "gravity_no_zeros.rdata", compress = "xz")

