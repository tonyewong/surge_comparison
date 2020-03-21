##==============================================================================
## calibration_driver.R
##
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================


rm(list=ls())

library(extRemes)
library(DEoptim)

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

##=============================================================================
## helper functions and some set up

# output file name
distr <- 'gev'
do_data_processing <- FALSE
today=Sys.Date(); today=format(today,format="%d%b%Y")
sim_id <- paste(distr,today, sep="-")
filename.optim <- paste("../output/optim_",sim_id,".rds", sep="")

if (do_data_processing) {
  source("process_gev.R")
  source("process_data.R") # do processing, then read the anticipated file names below
  data_calib <- readRDS(paste("../input_data/processeddata_gev_",today,".rds", sep=""))
  covariates <- readRDS(paste("../input_data/covariates_",today,".rds", sep=""))
} else {
  # read tide gauge data and covariates, fit previously
  data_calib <- readRDS("../input_data/processeddata_gev_21Mar2020.rds")
  covariates <- readRDS("../input_data/covariates_21Mar2020.rds")
}
site_names <- names(data_calib)
names_covariates <- colnames(covariates)[2:5]

# routine to trim the covariate forcing down to the length of the data set
source("trimmed_forcing.R")

# log-likelihood functions (prior and posterior not used here)
source("likelihood_gev.R")
#source("likelihood_gpd.R")

# set up parameters for all 8 candidate model structures for GEV and GPD
source("parameter_setup_gev.R")
#source("parameter_setup_gpd.R")
##=============================================================================



##=============================================================================
## MLE calibration for each tide gauge site

# settings for DEoptim (to find MCMC initial conditions)
NP.deoptim <- 100      # number of DE population members (at least 10*[# parameters])
niter.deoptim <- 100   # number of DE iterations
F.deoptim <- 0.8
CR.deoptim <- 0.9

optim_out <- vector("list", length(data_calib))
names(optim_out) <- names(data_calib)

for (dd in 1:length(data_calib)) {
  print(paste("Maximum likelihood estimation for site ",site_names[dd]," (",dd,"/",length(data_calib),")...", sep=""))
  optim_out[[dd]] <- vector("list", length(names_covariates))
  names(optim_out[[dd]]) <- names_covariates
  for (cc in names_covariates) {
    covar_forc <- covariates[,cc]
    time_forc <- covariates[,"year"]
    optim_out[[dd]][[cc]] <- vector("list", nmodel)
    for (mm in 1:nmodel) {
      if (mm > 1) {auxiliary <- trimmed_forcing(data_calib[[dd]][,"year"], time_forc, covar_forc)$forcing
      } else {auxiliary <- NULL}

      # initial parameter estimates
      out.deoptim <- DEoptim(neg_log_like_gev, lower=gev_models[[mm]]$bound_lower, upper=gev_models[[mm]]$bound_upper,
                             DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                             parnames=gev_models[[mm]]$parnames, data_calib=data_calib[[dd]][,"lsl_max"], auxiliary=auxiliary)
      optim_out[[dd]][[cc]][[mm]] <- out.deoptim
    }
  }
}
# save
print(paste('saving MCMC output as .rds file ', filename.optim, sep=''))
saveRDS(optim_out, file=filename.optim)
##=============================================================================



##==============================================================================
## End
##==============================================================================
