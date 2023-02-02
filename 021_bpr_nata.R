
# Title: Bayesian profile regression with NATA data
# Author: Kara E. McCormack
# Date Created: July 28, 2022
# Email: kem81@duke.edu

###################################################
### Preparations
###################################################
library(PReMiuM)
library(readr)
library(coda)
library(RColorBrewer)
library(dplyr)
library(readr)
library(here)
library(tidyverse)

# set.seed(1802)
seed <- 99488738

###################################################
### NATA data
###################################################

# read in data
inputs <- read_csv(here("02pol_dat_bpr.csv"))

# identify the SES variable names from data
nSweeps <- 500
nBurn <- 500
covNames <- colnames(inputs)
data <- inputs
output <- "nata_out"
nClusInit = c(1, 10, 30, 50) # number of initial clusters
plotBySweep = 25 # in final plots, how often (in terms of sweeps) to plot a point
# note: a ratio of about 1/10 seems to get good resolution if plotting
# a single plot.

##### A. Initial number of clusters vs. alpha ########
########## nClusInit = 1 ######
nClusInit = 1
# Run = 1
runInfoObj_nata1_1 <- profRegr(yModel = "Normal",
                              xModel = "Normal",
                              nSweeps = nSweeps,
                              nBurn = nBurn,
                              data = data,
                              output = output,
                              covNames = covNames,
                              nClusInit = nClusInit,
                              whichLabelSwitch="123",
                              run = TRUE,
                              excludeY = TRUE,
                              seed = 1)
# load information for diagnostics from output files
alpha = read.table(here("nata_out_alpha.txt"))
kappa1 = read.table(here("nata_out_kappa1.txt"))
logpost = read.table(here("nata_out_logPost.txt"))
nClusters = read.table(here("nata_out_nClusters.txt"))

res_test = cbind(sweep = 1:nSweeps,
                 alpha = alpha,
                 kappa = kappa1,
                 logPost = logpost[,1],
                 logLike = logpost[,2],
                 logPrior = logpost[,3],
                 nClusInit = rep(nClusInit, nSweeps),
                 nClusters = nClusters,
                 run = rep(1, length(nClusters)))
colnames(res_test)[2] <- "alpha"
colnames(res_test)[3] <- "kappa"
colnames(res_test)[8] <- "nClusters"
# save to results
r.nata_init1_1 <- res_test

# calculate clusters
disSimObj_nata1_1 <- calcDissimilarityMatrix(runInfoObj_nata1_1)
clusters_nata1_1 <- calcOptimalClustering(disSimObj_nata1_1,
                                      maxNClusters=NULL,
                                      useLS=F)
# nClusters_nata1_1 <- clusters_nata1_1$nClusters
# clustering_nata1_1 <- clusters_nata1_1$clustering

#----------------------
# nClusInit = 1
# Run = 2
runInfoObj_nata1_2 <- profRegr(yModel = "Normal",
                              xModel = "Normal",
                              nSweeps = nSweeps,
                              nBurn = nBurn,
                              data = data,
                              output = output,
                              covNames = covNames,
                              nClusInit = nClusInit,
                              whichLabelSwitch="123",
                              run = TRUE,
                              excludeY = TRUE,
                              seed = 2)
# load information for diagnostics from output files
alpha = read.table(here("nata_out_alpha.txt"))
kappa1 = read.table(here("nata_out_kappa1.txt"))
logpost = read.table(here("nata_out_logPost.txt"))
nClusters = read.table(here("nata_out_nClusters.txt"))

res_test = cbind(sweep = 1:nSweeps,
                 alpha = alpha,
                 kappa = kappa1,
                 logPost = logpost[,1],
                 logLike = logpost[,2],
                 logPrior = logpost[,3],
                 nClusInit = rep(nClusInit, nSweeps),
                 nClusters = nClusters,
                 run = rep(2, length(nClusters)))
colnames(res_test)[2] <- "alpha"
colnames(res_test)[3] <- "kappa"
colnames(res_test)[8] <- "nClusters"
# save to results
r.nata_init1_2 <- res_test

# calculate clusters
disSimObj_nata1_2 <- calcDissimilarityMatrix(runInfoObj_nata1_2)
clusters_nata1_2 <- calcOptimalClustering(disSimObj_nata1_2,
                                         maxNClusters=NULL,
                                         useLS=F)
nClusters_nata1_2 <- clusters_nata1_2$nClusters
clustering_nata1_2 <- clusters_nata1_2$clustering

#------------------------
# nClusInit = 1
# Run = 3
runInfoObj_nata1_3 <- profRegr(yModel = "Normal",
                              xModel = "Normal",
                              nSweeps = nSweeps,
                              nBurn = nBurn,
                              data = data,
                              output = output,
                              covNames = covNames,
                              nClusInit = nClusInit,
                              whichLabelSwitch="123",
                              run = TRUE,
                              excludeY = TRUE,
                              seed = 3)
# load information for diagnostics from output files
alpha = read.table(here("nata_out_alpha.txt"))
kappa1 = read.table(here("nata_out_kappa1.txt"))
logpost = read.table(here("nata_out_logPost.txt"))
nClusters = read.table(here("nata_out_nClusters.txt"))

res_test = cbind(sweep = 1:nSweeps,
                 alpha = alpha,
                 kappa = kappa1,
                 logPost = logpost[,1],
                 logLike = logpost[,2],
                 logPrior = logpost[,3],
                 nClusInit = rep(nClusInit, nSweeps),
                 nClusters = nClusters, 
                 run = rep(3, length(nClusters)))
colnames(res_test)[2] <- "alpha"
colnames(res_test)[3] <- "kappa"
colnames(res_test)[8] <- "nClusters"
# save to results
r.nata_init1_3 <- res_test

# calculate clusters
disSimObj_nata1_3 <- calcDissimilarityMatrix(runInfoObj_nata1_3)
clusters_nata1_3 <- calcOptimalClustering(disSimObj_nata1_3,
                                         maxNClusters=NULL,
                                         useLS=F)
nClusters_nata1_3 <- clusters_nata1_3$nClusters
clustering_nata1_3 <- clusters_nata1_3$clustering

#-------------------------------
####### nClusInit = 10 #######
nClusInit = 10
# Run = 1
runInfoObj_nata10_1 <- profRegr(yModel = "Normal",
                               xModel = "Normal",
                               nSweeps = nSweeps,
                               nBurn = nBurn,
                               data = data,
                               output = output,
                               covNames = covNames,
                               nClusInit = nClusInit,
                               whichLabelSwitch="123",
                               run = TRUE,
                               excludeY = TRUE,
                               seed = 43234)
# load information for diagnostics from output files
alpha = read.table(here("nata_out_alpha.txt"))
kappa1 = read.table(here("nata_out_kappa1.txt"))
logpost = read.table(here("nata_out_logPost.txt"))
nClusters = read.table(here("nata_out_nClusters.txt"))

res_test = cbind(sweep = 1:nSweeps,
                 alpha = alpha,
                 kappa = kappa1,
                 logPost = logpost[,1],
                 logLike = logpost[,2],
                 logPrior = logpost[,3],
                 nClusInit = rep(nClusInit, nSweeps),
                 nClusters = nClusters,
                 run = rep(1, length(nClusters)))
colnames(res_test)[2] <- "alpha"
colnames(res_test)[3] <- "kappa"
colnames(res_test)[8] <- "nClusters"
# save to results
r.nata_init10_1 <- res_test

# calculate clusters
disSimObj_nata10_1 <- calcDissimilarityMatrix(runInfoObj_nata10_1)
clusters_nata10_1 <- calcOptimalClustering(disSimObj_nata10_1,
                                         maxNClusters=NULL,
                                         useLS=F)
nClusters_nata10_1 <- clusters_nata10_1$nClusters
clustering_nata10_1 <- clusters_nata10_1$clustering

# calculate average risk profile
riskProfileObj10_1 <- calcAvgRiskAndProfile(clusters_nata10_1)

#----------------------
nClusInit = 10
# Run = 2
runInfoObj_nata10_2 <- profRegr(yModel = "Normal",
                               xModel = "Normal",
                               nSweeps = nSweeps,
                               nBurn = nBurn,
                               data = data,
                               output = output,
                               covNames = covNames,
                               nClusInit = nClusInit,
                               whichLabelSwitch="123",
                               run = TRUE,
                               excludeY = TRUE,
                               seed = 53425)
# load information for diagnostics from output files
alpha = read.table(here("nata_out_alpha.txt"))
kappa1 = read.table(here("nata_out_kappa1.txt"))
logpost = read.table(here("nata_out_logPost.txt"))
nClusters = read.table(here("nata_out_nClusters.txt"))

res_test = cbind(sweep = 1:nSweeps,
                 alpha = alpha,
                 kappa = kappa1,
                 logPost = logpost[,1],
                 logLike = logpost[,2],
                 logPrior = logpost[,3],
                 nClusInit = rep(nClusInit, nSweeps),
                 nClusters = nClusters,
                 run = rep(2, length(nClusters)))
colnames(res_test)[2] <- "alpha"
colnames(res_test)[3] <- "kappa"
colnames(res_test)[8] <- "nClusters"
# save to results
results_init10_2 <- res_test

# calculate clusters
disSimObj_nata10_2 <- calcDissimilarityMatrix(runInfoObj_nata10_2)
clusters_nata10_2 <- calcOptimalClustering(disSimObj_nata10_2,
                                          maxNClusters=NULL,
                                          useLS=F)
nClusters_nata10_2 <- clusters_nata10_2$nClusters
clustering_nata10_2 <- clusters_nata10_2$clustering

# calculate average risk profile
riskProfileObj10_2 <- calcAvgRiskAndProfile(clusters_nata10_2)

#-----------------------

nClusInit = 10
# Run = 3
runInfoObj_nata10_3 <- profRegr(yModel = "Normal",
                               xModel = "Normal",
                               nSweeps = nSweeps,
                               nBurn = nBurn,
                               data = data,
                               output = output,
                               covNames = covNames,
                               nClusInit = nClusInit,
                               whichLabelSwitch="123",
                               run = TRUE,
                               excludeY = TRUE,
                               seed = 64445)
# load information for diagnostics from output files
alpha = read.table(here("nata_out_alpha.txt"))
kappa1 = read.table(here("nata_out_kappa1.txt"))
logpost = read.table(here("nata_out_logPost.txt"))
nClusters = read.table(here("nata_out_nClusters.txt"))

res_test = cbind(sweep = 1:nSweeps,
                 alpha = alpha,
                 kappa = kappa1,
                 logPost = logpost[,1],
                 logLike = logpost[,2],
                 logPrior = logpost[,3],
                 nClusInit = rep(nClusInit, nSweeps),
                 nClusters = nClusters, 
                 run = rep(3, length(nClusters)))
colnames(res_test)[2] <- "alpha"
colnames(res_test)[3] <- "kappa"
colnames(res_test)[8] <- "nClusters"
# save to results
r.nata_init10_3 <- res_test    

# calculate clusters
disSimObj_nata10_3 <- calcDissimilarityMatrix(runInfoObj_nata10_3)
clusters_nata10_3 <- calcOptimalClustering(disSimObj_nata10_3,
                                          maxNClusters=NULL,
                                          useLS=F)
nClusters_nata10_3 <- clusters_nata10_3$nClusters
clustering_nata10_3 <- clusters_nata10_3$clustering

# calculate average risk profile
riskProfileObj10_3 <- calcAvgRiskAndProfile(clusters_acs10_3)


####### nClusInit = 30    #####
nClusInit = 30
# Run = 1
runInfoObj_nata30_1 <- profRegr(yModel = "Normal",
                               xModel = "Normal",
                               nSweeps = nSweeps,
                               nBurn = nBurn,
                               data = data,
                               output = output,
                               covNames = covNames,
                               nClusInit = nClusInit,
                               whichLabelSwitch="123",
                               run = TRUE,
                               excludeY = TRUE,
                               seed = 7667)
# load information for diagnostics from output files
alpha = read.table(here("nata_out_alpha.txt"))
kappa1 = read.table(here("nata_out_kappa1.txt"))
logpost = read.table(here("nata_out_logPost.txt"))
nClusters = read.table(here("nata_out_nClusters.txt"))

res_test = cbind(sweep = 1:nSweeps,
                 alpha = alpha,
                 kappa = kappa1,
                 logPost = logpost[,1],
                 logLike = logpost[,2],
                 logPrior = logpost[,3],
                 nClusInit = rep(nClusInit, nSweeps),
                 nClusters = nClusters,
                 run = rep(1, length(nClusters)))
colnames(res_test)[2] <- "alpha"
colnames(res_test)[3] <- "kappa"
colnames(res_test)[8] <- "nClusters"
# save to results
results_init30_1 <- res_test

# calculate clusters
disSimObj_nata30_1 <- calcDissimilarityMatrix(runInfoObj_nata30_1)
clusters_nata30_1 <- calcOptimalClustering(disSimObj_nata30_1,
                                          maxNClusters=NULL,
                                          useLS=F)
nClusters_nata30_1 <- clusters_nata30_1$nClusters
clustering_nata30_1 <- clusters_nata30_1$clustering

# calculate average risk profile
riskProfileObj30_1 <- calcAvgRiskAndProfile(clusters_nata30_1)

#-----------------

nClusInit = 30
# Run = 2
runInfoObj_nata30_2 <- profRegr(yModel = "Normal",
                               xModel = "Normal",
                               nSweeps = nSweeps,
                               nBurn = nBurn,
                               data = data,
                               output = output,
                               covNames = covNames,
                               nClusInit = nClusInit,
                               whichLabelSwitch="123",
                               run = TRUE,
                               excludeY = TRUE,
                               seed = 835656)
# load information for diagnostics from output files
alpha = read.table(here("nata_out_alpha.txt"))
kappa1 = read.table(here("nata_out_kappa1.txt"))
logpost = read.table(here("nata_out_logPost.txt"))
nClusters = read.table(here("nata_out_nClusters.txt"))

res_test = cbind(sweep = 1:nSweeps,
                 alpha = alpha,
                 kappa = kappa1,
                 logPost = logpost[,1],
                 logLike = logpost[,2],
                 logPrior = logpost[,3],
                 nClusInit = rep(nClusInit, nSweeps),
                 nClusters = nClusters,
                 run = rep(2, length(nClusters)))
colnames(res_test)[2] <- "alpha"
colnames(res_test)[3] <- "kappa"
colnames(res_test)[8] <- "nClusters"
# save to results
r.nata_init30_2 <- res_test

# calculate clusters
disSimObj_nata30_2 <- calcDissimilarityMatrix(runInfoObj_nata30_2)
clusters_nata30_2 <- calcOptimalClustering(disSimObj_nata30_2,
                                          maxNClusters=NULL,
                                          useLS=F)
nClusters_nata30_2 <- clusters_nata30_2$nClusters
clustering_nata30_2 <- clusters_nata30_2$clustering

# calculate average risk profile
riskProfileObj30_2 <- calcAvgRiskAndProfile(clusters_nata30_2)

#-----------------------------
# Run = 3
runInfoObj_nata30_3 <- profRegr(yModel = "Normal",
                               xModel = "Normal",
                               nSweeps = nSweeps,
                               nBurn = nBurn,
                               data = data,
                               output = output,
                               covNames = covNames,
                               nClusInit = nClusInit,
                               whichLabelSwitch="123",
                               run = TRUE,
                               excludeY = TRUE,
                               seed = 9344)
# load information for diagnostics from output files
alpha = read.table(here("nata_out_alpha.txt"))
kappa1 = read.table(here("nata_out_kappa1.txt"))
logpost = read.table(here("nata_out_logPost.txt"))
nClusters = read.table(here("nata_out_nClusters.txt"))

res_test = cbind(sweep = 1:nSweeps,
                 alpha = alpha,
                 kappa = kappa1,
                 logPost = logpost[,1],
                 logLike = logpost[,2],
                 logPrior = logpost[,3],
                 nClusInit = rep(nClusInit, nSweeps),
                 nClusters = nClusters, 
                 run = rep(3, length(nClusters)))
colnames(res_test)[2] <- "alpha"
colnames(res_test)[3] <- "kappa"
colnames(res_test)[8] <- "nClusters"
# save to results
r.nata_init30_3 <- res_test    

# calculate clusters
disSimObj_nata30_3 <- calcDissimilarityMatrix(runInfoObj_nata30_3)
clusters_nata30_3 <- calcOptimalClustering(disSimObj_nata30_3,
                                          maxNClusters=NULL,
                                          useLS=F)
nClusters_nata30_3 <- clusters_nata30_3$nClusters
clustering_nata30_3 <- clusters_nata30_3$clustering

# calculate average risk profile
riskProfileObj30_3 <- calcAvgRiskAndProfile(clusters_nata30_3)

#######----- nClusInit = 50    #####    
nClusInit = 50
# Run = 1
runInfoObj_nata50_1 <- profRegr(yModel = "Normal",
                               xModel = "Normal",
                               nSweeps = nSweeps,
                               nBurn = nBurn,
                               data = data,
                               output = output,
                               covNames = covNames,
                               nClusInit = nClusInit,
                               whichLabelSwitch="123",
                               run = TRUE,
                               excludeY = TRUE,
                               seed = 1066)
# load information for diagnostics from output files
alpha = read.table(here("nata_out_alpha.txt"))
kappa1 = read.table(here("nata_out_kappa1.txt"))
logpost = read.table(here("nata_out_logPost.txt"))
nClusters = read.table(here("nata_out_nClusters.txt"))

res_test = cbind(sweep = 1:nSweeps,
                 alpha = alpha,
                 kappa = kappa1,
                 logPost = logpost[,1],
                 logLike = logpost[,2],
                 logPrior = logpost[,3],
                 nClusInit = rep(nClusInit, nSweeps),
                 nClusters = nClusters,
                 run = rep(1, length(nClusters)))
colnames(res_test)[2] <- "alpha"
colnames(res_test)[3] <- "kappa"
colnames(res_test)[8] <- "nClusters"
# save to results
r.nata_init50_1 <- res_test

# calculate clusters
disSimObj_nata50_1 <- calcDissimilarityMatrix(runInfoObj_nata50_1)
clusters_nata50_1 <- calcOptimalClustering(disSimObj_nata50_1,
                                          maxNClusters=NULL,
                                          useLS=F)
nClusters_nata50_1 <- clusters_nata50_1$nClusters
clustering_nata50_1 <- clusters_nata50_1$clustering

# calculate average risk profile
riskProfileObj50_1 <- calcAvgRiskAndProfile(clusters_nata50_1)

#-------------------------------

# Run = 2
runInfoObj_nata50_2 <- profRegr(yModel = "Normal",
                               xModel = "Normal",
                               nSweeps = nSweeps,
                               nBurn = nBurn,
                               data = data,
                               output = output,
                               covNames = covNames,
                               nClusInit = nClusInit,
                               whichLabelSwitch="123",
                               run = TRUE,
                               excludeY = TRUE,
                               seed = 11544)
# load information for diagnostics from output files
alpha = read.table(here("nata_out_alpha.txt"))
kappa1 = read.table(here("nata_out_kappa1.txt"))
logpost = read.table(here("nata_out_logPost.txt"))
nClusters = read.table(here("nata_out_nClusters.txt"))

res_test = cbind(sweep = 1:nSweeps,
                 alpha = alpha,
                 kappa = kappa1,
                 logPost = logpost[,1],
                 logLike = logpost[,2],
                 logPrior = logpost[,3],
                 nClusInit = rep(nClusInit, nSweeps),
                 nClusters = nClusters,
                 run = rep(2, length(nClusters)))
colnames(res_test)[2] <- "alpha"
colnames(res_test)[3] <- "kappa"
colnames(res_test)[8] <- "nClusters"
# save to results
r.nata_init50_2 <- res_test

# calculate clusters
disSimObj_nata50_2 <- calcDissimilarityMatrix(runInfoObj_nata50_2)
clusters_nata50_2 <- calcOptimalClustering(disSimObj_nata50_2,
                                          maxNClusters=NULL,
                                          useLS=F)
nClusters_nata50_2 <- clusters_nata50_2$nClusters
clustering_nata50_2 <- clusters_nata50_2$clustering

# calculate average risk profile
riskProfileObjnata50_2 <- calcAvgRiskAndProfile(clusters_nata50_2)

#----------------------
# Run = 3
runInfoObj_nata50_3 <- profRegr(yModel = "Normal",
                               xModel = "Normal",
                               nSweeps = nSweeps,
                               nBurn = nBurn,
                               data = data,
                               output = output,
                               covNames = covNames,
                               nClusInit = nClusInit,
                               whichLabelSwitch="123",
                               run = TRUE,
                               excludeY = TRUE,
                               seed = 1265656)
# load information for diagnostics from output files
alpha = read.table(here("nata_out_alpha.txt"))
kappa1 = read.table(here("nata_out_kappa1.txt"))
logpost = read.table(here("nata_out_logPost.txt"))
nClusters = read.table(here("nata_out_nClusters.txt"))

res_test = cbind(sweep = 1:nSweeps,
                 alpha = alpha,
                 kappa = kappa1,
                 logPost = logpost[,1],
                 logLike = logpost[,2],
                 logPrior = logpost[,3],
                 nClusInit = rep(nClusInit, nSweeps),
                 nClusters = nClusters, 
                 run = rep(3, length(nClusters)))
colnames(res_test)[2] <- "alpha"
colnames(res_test)[3] <- "kappa"
colnames(res_test)[8] <- "nClusters"
# save to results
r.nata_init50_3 <- res_test      

# calculate clusters
disSimObj_nata50_3 <- calcDissimilarityMatrix(runInfoObj_nata50_3)
clusters_nata50_3 <- calcOptimalClustering(disSimObj_nata50_3,
                                          maxNClusters=NULL,
                                          useLS=F)
nClusters_nata50_3 <- clusters_nata50_3$nClusters
clustering_nata50_3 <- clusters_nata50_3$clustering

# calculate average risk profile
riskProfileObjnata50_3 <- calcAvgRiskAndProfile(clusters_nata50_3)

#-----------------------------

###### plot diagnostic A #######
res_acs <- rbind(r.nata_init1_1, r.nata_init1_2, r.nata_init1_3, 
                 r.nata_init10_1, r.nata_init10_2, r.nata_init10_3, 
                 r.nata_init30_1, r.nata_init30_2, r.nata_init30_3, 
                 r.nata_init50_1, r.nata_init50_2, r.nata_init50_3)

pA_nata <- res_nata %>%
  ggplot(aes(x = as.factor(nClusInit), y = alpha, fill = as.factor(run))) +
  geom_boxplot() +
  ggtitle("NATA Data")

###### B. Initial number of clusters vs. posterior number of clusters. ######

# just need to run it once more with nClusInit = 5
nClusInit <- 5
runInfoObj_nata5_1 <- profRegr(yModel = "Normal",
                              xModel = "Normal",
                              nSweeps = nSweeps,
                              nBurn = nBurn,
                              data = data,
                              output = output,
                              covNames = covNames,
                              nClusInit = nClusInit,
                              whichLabelSwitch="123",
                              run = TRUE,
                              excludeY = TRUE,
                              seed = 13690) # change seed
# load information for diagnostics from output files
alpha = read.table(here("nata_out_alpha.txt"))
kappa1 = read.table(here("nata_out_kappa1.txt"))
logpost = read.table(here("nata_out_logPost.txt"))
nClusters = read.table(here("nata_out_nClusters.txt"))

res_test = cbind(sweep = 1:nSweeps,
                 alpha = alpha,
                 kappa = kappa1,
                 logPost = logpost[,1],
                 logLike = logpost[,2],
                 logPrior = logpost[,3],
                 nClusInit = rep(nClusInit, nSweeps),
                 nClusters = nClusters, 
                 run = rep(1, length(nClusters))) # change to run 1
colnames(res_test)[2] <- "alpha"
colnames(res_test)[3] <- "kappa"
colnames(res_test)[8] <- "nClusters"

# save to results
r.nata_init5_1 <- res_test   

# concatenate first run from each initial # clusters
resB_nata <- rbind(r.nata_init1_1,
                  r.nata_init5_1,
                  r.nata_init10_1, 
                  r.nata_init30_1, 
                  r.nata_init50_1)

pB_nata <- resB_nata %>%
  ggplot(aes(x = as.factor(nClusInit), y = nClusters)) +
  geom_boxplot() +
  ggtitle("NATA Data")

# save .RData
save.image(here("resultsA_nata.RData"))

print("hello world")


