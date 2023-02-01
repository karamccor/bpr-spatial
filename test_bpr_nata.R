
# read in data
inputs <- read_csv(here("data/02pol_dat_bpr.csv"))

# global parameters

covNames <- colnames(inputs) # identify variable names from columns of data
data <- inputs
output <- "nata_out" # overrwrite output files each time, will help with storage

runInfoObj_nata1 <- profRegr(yModel = "Normal", 
                             xModel = "Normal", 
                             nSweeps = 500, 
                             nBurn = 500,
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