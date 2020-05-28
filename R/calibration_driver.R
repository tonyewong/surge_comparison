##==============================================================================
## calibration_driver.R
##
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================


rm(list=ls())

library(extRemes)
library(DEoptim)
library(date)

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

##==============================================================================
## Helpers and set up

calib_post <- TRUE
do_data_processing <- FALSE
do_fit_priors <- FALSE
min_years <- 15 # minimum number of years for using a tide gauge site
today=Sys.Date(); today=format(today,format="%d%b%Y")

# helper functions
source("trimmed_forcing.R")
source("decluster_timeseries.R")
source("process_gev.R")
source("process_gpd.R")

# settings for DEoptim (to minimize the negative log-likelihood)
NP.deoptim <- 100      # number of DE population members (at least 10*[# parameters])
niter.deoptim <- 100   # number of DE iterations
F.deoptim <- 0.8
CR.deoptim <- 0.9
##==============================================================================



##=============================================================================
## data processing for GEV, GPD and covariates

if (do_data_processing) {
  source("process_data.R") # do processing, then read the anticipated file names below
  data_gev <- readRDS(paste("../input_data/processeddata_gev_",today,".rds", sep=""))
  data_gpd <- readRDS(paste("../input_data/processeddata_gpd_",today,".rds", sep=""))
  covariates <- readRDS(paste("../input_data/covariates_",today,".rds", sep=""))
} else {
  # read tide gauge data and covariates, fit previously
  data_gev <- readRDS("../input_data/processeddata_gev_22Mar2020.rds")
  data_gpd <- readRDS("../input_data/processeddata_gpd_23Mar2020.rds")
  covariates <- readRDS("../input_data/covariates_22Mar2020.rds")
}
site_names <- names(data_gev)
names_covariates <- colnames(covariates)[2:ncol(covariates)]
##==============================================================================



##=============================================================================
## priors

if (do_fit_priors) {
  source("fit_priors.R") # do processing, then read the anticipated file names below
} else {
  # read priors, fit previously
  priors <- vector("list", 2); names(priors) <- c("gev","gpd")
  priors$gev <- readRDS("../input_data/surge_priors_normalgamma_gev_15May2020.rds")
  priors$gpd <- readRDS("../input_data/surge_priors_normalgamma_gpd_15May2020.rds")
}
##==============================================================================



##==============================================================================
## MLE or MAP calibration for GEV parameters...
## ... for each model, for each covariate, for each tide gauge site

distr <- 'gev'
if (calib_post) {sim_id <- paste(distr,"post",today, sep="_")
} else {sim_id <- paste(distr,"like",today, sep="_")}
filename.optim <- paste("../output/optim_",sim_id,".rds", sep="")

# log-likelihood functions (prior and posterior not used here)
source("likelihood_gev.R")

# set up parameters for all 8 candidate model structures
source("parameter_setup_gev.R")

optim_out <- vector("list", length(data_gev))
names(optim_out) <- names(data_gev)

for (dd in 1:length(data_gev)) {
  print(paste("Parameter estimation for GEV models for site ",site_names[dd]," (",dd,"/",length(data_gev),")...", sep=""))
  optim_out[[dd]] <- vector("list", length(names_covariates))
  names(optim_out[[dd]]) <- names_covariates
  for (cc in names_covariates) {
    covar_forc <- covariates[,cc]
    time_forc <- covariates[,"year"]
    optim_out[[dd]][[cc]] <- vector("list", nmodel)
    for (mm in 1:nmodel) {
      if (mm > 1) {auxiliary <- trimmed_forcing(data_gev[[dd]][,"year"], time_forc, covar_forc)$forcing
      } else {auxiliary <- NULL}
      if (calib_post) {
        out.deoptim <- DEoptim(neg_log_post_gev, lower=gev_models[[mm]]$bound_lower, upper=gev_models[[mm]]$bound_upper,
                               DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                               parnames=gev_models[[mm]]$parnames, data_calib=data_gev[[dd]][,"lsl_max"],
                               priors=priors$gev[[cc]][[mm]], auxiliary=auxiliary)
      } else {
        out.deoptim <- DEoptim(neg_log_like_gev, lower=gev_models[[mm]]$bound_lower, upper=gev_models[[mm]]$bound_upper,
                               DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                               parnames=gev_models[[mm]]$parnames, data_calib=data_gev[[dd]][,"lsl_max"], auxiliary=auxiliary)
      }
      optim_out[[dd]][[cc]][[mm]] <- out.deoptim
    }
  }
  # takes a while, so save after each one
  saveRDS(optim_out, file=filename.optim)
}
print(paste('Saved DEoptim output as .rds file ', filename.optim, sep=''))
##==============================================================================



##==============================================================================
## MLE or MAP calibration for GPD parameters...
## ... for each model, for each covariate, for each tide gauge site

distr <- 'gpd'
if (calib_post) {sim_id <- paste(distr,"post",today, sep="_")
} else {sim_id <- paste(distr,"like",today, sep="_")}
filename.optim <- paste("../output/optim_",sim_id,".rds", sep="")

# log-likelihood functions (prior and posterior not used here)
source("likelihood_gpd.R")

# set up parameters for all 8 candidate model structures
source("parameter_setup_gpd.R")

optim_out <- vector("list", length(data_gpd))
names(optim_out) <- names(data_gpd)

for (dd in 1:length(data_gpd)) {
  print(paste("Parameter estimation for GPD models for site ",site_names[dd]," (",dd,"/",length(data_gpd),")...", sep=""))
  optim_out[[dd]] <- vector("list", length(names_covariates))
  names(optim_out[[dd]]) <- names_covariates
  for (cc in names_covariates) {
    covar_forc <- covariates[,cc]
    time_forc <- covariates[,"year"]
    optim_out[[dd]][[cc]] <- vector("list", nmodel)
    for (mm in 1:nmodel) {
      if (mm > 1) {auxiliary <- trimmed_forcing(data_gpd[[dd]]$year, time_forc, covar_forc)$forcing
      } else {auxiliary <- NULL}
      if (calib_post) {
        out.deoptim <- DEoptim(neg_log_post_ppgpd, lower=gpd_models[[mm]]$bound_lower, upper=gpd_models[[mm]]$bound_upper,
                               DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                               parnames=gpd_models[[mm]]$parnames, data_calib=data_gpd[[dd]],
                               priors=priors$gpd[[cc]], model=mm, auxiliary=auxiliary)
      } else {
        out.deoptim <- DEoptim(neg_log_like_ppgpd, lower=gpd_models[[mm]]$bound_lower, upper=gpd_models[[mm]]$bound_upper,
                               DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                               parnames=gpd_models[[mm]]$parnames, data_calib=data_gpd[[dd]], auxiliary=auxiliary)
      }
      optim_out[[dd]][[cc]][[mm]] <- out.deoptim
    }
  }
  # takes a while, so save after each one
  saveRDS(optim_out, file=filename.optim)
}
print(paste('Saved DEoptim output as .rds file ', filename.optim, sep=''))
##==============================================================================



##==============================================================================
## End
##==============================================================================
