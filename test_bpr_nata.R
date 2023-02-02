
# load libraries
library(readr)
library(tidyverse)
library(PReMiuM)

# read in data
# select three pollutants for now
inputs <- read_csv("~/school/aim1-spatial/data/02pol_dat_bpr.csv") %>%
  # select(`1,3-BUTADIENE`, `ACETALDEHYDE`, `BENZENE`, 
  #        `DIESEL PM`, `ETHYLBENZENE`, `FORMALDEHYDE`, `HEXANE`, `LEAD COMPOUNDS`, 
  #        `MANGANESE COMPOUNDS`, `MERCURY COMPOUNDS`) %>%
  rename(butatiene = "1,3-BUTADIENE",
         acetaldehyde = "ACETALDEHYDE",
         benzene = "BENZENE", 
         diesel_pm = "DIESEL PM",
         ethylbenzene = "ETHYLBENZENE",
         formaldehyde = "FORMALDEHYDE", 
         hexane = "HEXANE", 
         lead = "LEAD COMPOUNDS", 
         manganese = "MANGANESE COMPOUNDS", 
         mercury = "MERCURY COMPOUNDS", 
         methanol = "METHANOL", 
         nickel = "NICKEL COMPOUNDS", 
         toluene = "TOLUENE", 
         xylenes = "XYLENES (MIXED ISOMERS)")

hist(inputs$`1,3-BUTADIENE`)
hist(inputs$BENZENE)
hist(inputs$`DIESEL PM`)
hist(inputs$METHANOL)
hist(inputs$HEXANE)
hist(inputs$`LEAD COMPOUNDS`)

# global parameters

covNames <- colnames(inputs) # identify variable names from columns of data
data <- inputs
output <- "nata_out" # overrwrite output files each time, will help with storage

runInfoObj_nata1 <- profRegr(yModel = "Normal", 
                             xModel = "Normal", 
                             nSweeps = 100, 
                             nBurn = 100,
                             data = data, 
                             output = output, 
                             covNames = covNames,  
                             nClusInit = 5,
                             whichLabelSwitch="123", 
                             run = TRUE, excludeY = TRUE, seed = 123)

margModelPosterior(runInfoObj_acs1)

globalParsTrace(runInfoObj_acs1, parameters = "nClusters", 
                plotBurnIn = F, whichBeta = 1)

globalParsTrace(runInfoObj_acs1, parameters = "alpha", 
                plotBurnIn = F, whichBeta = 1)