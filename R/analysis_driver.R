##==============================================================================
## analysis_driver.R
##
## To be run after calibration_driver.R
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

rm(list=ls())

library(extRemes)
library(Bolstad)

if(Sys.info()['user']=='tony') {
  # Tony's local machine (if you aren't me, you almost certainly need to change this...)
  machine <- 'local'
  setwd('/Users/tony/codes/surge_comparison/R')
} else if(Sys.info()['user']=='aewsma') {
  machine <- 'office'
  setwd('/Users/aewsma/codes/surge_comparison/R')
} else {
  # assume on another cluster of some kind...
  machine <- 'remote'
  setwd('~/work/codes/surge_comparison/R')
}

## Set up. Make sure consistent with the
distr <- "gev"
datestamp <- "22Mar2020"
sim_id <- paste(distr,datestamp, sep="-")
filename_optim <- paste("../output/optim_",sim_id,".rds", sep="")
filename_covariates <- "../input_data/covariates_22Mar2020.rds"
filename_calibdata  <- "../input_data/processeddata_gev_22Mar2020.rds"

## do some analyzing
source("best_models.R") # yields parameters[[model]][site,parameter]
source("make_projections.R") # yields return level estimates

## number of years of data for each site
num_years <- rep(NA, nsite); names(num_years) <- site_names
for (dd in 1:nsite) {num_years[dd] <- nrow(data_calib[[dd]])}


## get to the good stuff!



##==============================================================================
## Compare model choices based
## on number of parameters
##============================

todo -- rework this into a single figure with 2 columns, and need to address covariate uncertainty

todo -- for covariate uncertainty , create a second similar figure that shows how choice of covariate is affected

todo -- then have a second version where instead of NLL/AIC/BIC, which evaluates number of parameters effect, show differences based on amount of available data.
     -- maybe bin into 3 roughly equally sized bins somehow? 15-45, 45-70, 70+?
     > length(which(num_years <= 45))
     [1] 12
     > length(which((num_years > 45) & (num_years <= 70)))
     [1] 12
     > length(which(num_years > 70))
     [1] 12
idx_short <- which(num_years <= 45)
idx_medium <- which((num_years > 45) & (num_years <= 70))
idx_long <- which(num_years > 70)

##  comparison of model choice based on penalization of over-fitting
best_models <- vector("list", 3); names(best_models) <- c("nll","aic","bic")
model_choices <- c(1,2,3,4,5,6,7,8)
par(mfrow=c(3,1))

metric <- "nll"
best_models[[metric]] <- rep(NA, nsite)
names(best_models[[metric]]) <- site_names
for (dd in site_names) {
    mat <- nll[dd,,model_choices]
    idx_min <- which.min(mat)
    idx_row <- idx_min%%nrow(mat)
    if(idx_row==0) {idx_row <- nrow(mat)}
    idx_col <- which.min(as.vector(mat[idx_row,]))
    best_models[[metric]][dd] <- idx_col
}
hist(best_models[[metric]], breaks=seq(from=0.5, to=8.5, by=1), xlim=c(0.5,8.5), main=metric)

metric <- "aic"
best_models[[metric]] <- rep(NA, nsite)
names(best_models[[metric]]) <- site_names
for (dd in site_names) {
    mat <- aic[dd,,model_choices]
    idx_min <- which.min(mat)
    idx_row <- idx_min%%nrow(mat)
    if(idx_row==0) {idx_row <- nrow(mat)}
    idx_col <- which.min(as.vector(mat[idx_row,]))
    best_models[[metric]][dd] <- idx_col
}
hist(best_models[[metric]], breaks=seq(from=0.5, to=8.5, by=1), xlim=c(0.5,8.5), main=metric)

metric <- "bic"
best_models[[metric]] <- rep(NA, nsite)
names(best_models[[metric]]) <- site_names
for (dd in site_names) {
    mat <- bic[dd,,model_choices]
    idx_min <- which.min(mat)
    idx_row <- idx_min%%nrow(mat)
    if(idx_row==0) {idx_row <- nrow(mat)}
    idx_col <- which.min(as.vector(mat[idx_row,]))
    best_models[[metric]][dd] <- idx_col
}
hist(best_models[[metric]], breaks=seq(from=0.5, to=8.5, by=1), xlim=c(0.5,8.5), main=metric)

##  same, but restrict to only stationary, mu nonstat, or sigma nonstat
best_models <- vector("list", 3); names(best_models) <- c("nll","aic","bic")
model_choices <- c(1,2,3)
par(mfrow=c(3,1))

metric <- "nll"
best_models[[metric]] <- rep(NA, nsite)
names(best_models[[metric]]) <- site_names
for (dd in site_names) {
    mat <- nll[dd,,model_choices]
    idx_min <- which.min(mat)
    idx_row <- idx_min%%nrow(mat)
    if(idx_row==0) {idx_row <- nrow(mat)}
    idx_col <- which.min(as.vector(mat[idx_row,]))
    best_models[[metric]][dd] <- idx_col
}
hist(best_models[[metric]], breaks=seq(from=0.5, to=8.5, by=1), xlim=c(0.5,8.5), main=metric)

metric <- "aic"
best_models[[metric]] <- rep(NA, nsite)
names(best_models[[metric]]) <- site_names
for (dd in site_names) {
    mat <- aic[dd,,model_choices]
    idx_min <- which.min(mat)
    idx_row <- idx_min%%nrow(mat)
    if(idx_row==0) {idx_row <- nrow(mat)}
    idx_col <- which.min(as.vector(mat[idx_row,]))
    best_models[[metric]][dd] <- idx_col
}
hist(best_models[[metric]], breaks=seq(from=0.5, to=8.5, by=1), xlim=c(0.5,8.5), main=metric)

metric <- "bic"
best_models[[metric]] <- rep(NA, nsite)
names(best_models[[metric]]) <- site_names
for (dd in site_names) {
    mat <- bic[dd,,model_choices]
    idx_min <- which.min(mat)
    idx_row <- idx_min%%nrow(mat)
    if(idx_row==0) {idx_row <- nrow(mat)}
    idx_col <- which.min(as.vector(mat[idx_row,]))
    best_models[[metric]][dd] <- idx_col
}
hist(best_models[[metric]], breaks=seq(from=0.5, to=8.5, by=1), xlim=c(0.5,8.5), main=metric)
##==============================================================================



##==============================================================================
## Compare model choices based
## on number of data available
##============================

todo

x axis -- number of years of block max data
y axis -- model choice (# parameters, or just number? so we can tell which structures are favored)
 -- use AIC as a happy medium?

##==============================================================================



##==============================================================================
## Compare return levels based
## on model choice
##============================

todo

##==============================================================================



##==============================================================================
##

##==============================================================================



##==============================================================================
## End
##==============================================================================
