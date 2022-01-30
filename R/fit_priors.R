##===============================================================================
## fit_priors.R
##
## Fit prior distributions for eight candidate GEV and PP/GPD-based model structures,
## times four candidate covariates, using all UHSLC data (research quality
## versions) at least 90 years, plus Sewells Point (Norfolk).
## Makes for 29 stations total.
##
## 1. get tide gauge data objects for all UHSLC database sites with > 90 years
## 2. get ... for Norfolk, VA, USA (not in database)
## 3. calculate maximum likelihood pp/gpd parameters, for each of the candidate
##    model structures, for each of the sites
## 4. fit normal or gamma prior distributions to these parameter sets, for each
##    model parameter within each of the candidate model structures.
## 5. write this priors object to a file (rds) and save progress to revisit later
##    (rdata)
##
## Updated 11 Dec 2017 // revised processing // tony wong
## Updated 18 Aug 2018 // revised for covariates // tony wong
## Updated 14 May 2020 // revised for comparison of extreme value models // tony wong
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.
##==============================================================================


rm(list=ls())

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

pot.threshold <- 0.99   # POT threshold (percentile)
dt.decluster <- 3       # declustering time-scale (days)

.NP.deoptim <- 100      # number of DE population members (at least 10*[# parameters])
.niter.deoptim <- 100   # number of DE iterations
output.dir <- '../input_data/'

filename.many <- '../input_data/tidegauge_processed_manystations_decl3-pot99-annual_10Dec2017.rds'
priors <- vector("list", 2); names(priors) <- c("gev","gpd")

#
#===============================================================================
# relevant libraries - do 'install.packages([library name])' if you do not have
# one yet
#===============================================================================
#

library(extRemes)
library(DEoptim)

#
#===============================================================================
# read and process data for forcing (auxiliary covariate for nonstationary
# parameters).
#===============================================================================
#

print('reading forcing data...')

# Load the covariate names and values
source('get_timeseries_covariates.R')

print('...done.')

#
#===============================================================================
# read and process data for tide gauge stations (or read previous processing)
#===============================================================================
#

print('reading processed data from tide gauge stations...')

data_many <- readRDS(filename.many)

# remove RQ#240 (Ferdinanda Beach, USA) because it's missing about 60 years
# from 1924-1984
data_many$rqh0240 <- NULL

# round them all up as one big data set
data_all <- data_many
# any others to add?

print('...done.')

#
#===============================================================================
# set up model parameters
#===============================================================================
#

print('setting up GEV and PP-GPD model parameters for DE optimization...')

source('parameter_setup_gev.R')
source('parameter_setup_gpd.R')

print('...done.')

#
#===============================================================================
# parameters for DE optim (for maximum likelihood/minimum negative likelihoood)
#===============================================================================
#

NP.deoptim <- .NP.deoptim
niter.deoptim <- .niter.deoptim
F.deoptim <- 0.8
CR.deoptim <- 0.9



#===============================================================================
#===============================================================================
# GPD
#===============================================================================
#===============================================================================



#
#===============================================================================
# fit MLE PP-GPD model parameters for each candidate model at each tide gauge
#===============================================================================
#

print('starting DE optimization for MLE PP-GPD parameters for all stations in set...')

# need the likelihood function
source('likelihood_gpd.R')

deoptim.all <- vector('list', length(names_covariates))
names(deoptim.all) <- names_covariates
for (cc in names_covariates) {
  deoptim.all[[cc]] <- vector('list', length(gpd_models))
  for (i in 1:length(gpd_models)) {
    deoptim.all[[cc]][[i]] <- mat.or.vec(length(data_all), length(gpd_models[[i]]$parnames))
    rownames(deoptim.all[[cc]][[i]]) <- names(data_all)
    colnames(deoptim.all[[cc]][[i]]) <- gpd_models[[i]]$parnames
  }
}

for (dd in 1:length(data_all)) {

  print(paste('starting to calculate MLE PP-GPD parameters for tide gauge data set ',dd,' / ',length(data_all),sep=''))
  tbeg0 <- proc.time()
  data_all[[dd]]$deoptim <- vector('list', nmodel)

  for (gpd.type in 1:length(gpd_models)) {

    print(paste('  - starting DE optimization for model',gpd.type,'...'))
    tbeg <- proc.time()

    irem_aux <- NULL

    # if tide gauge record starts before auxiliary forcing, clip it
    if(data_all[[dd]]$gpd$year[1] < covariates[1,'year']) {
      irem <- which(data_all[[dd]]$gpd$year < covariates[,'year'][1])
      data_all[[dd]]$gev_year$year <- data_all[[dd]]$gev_year$year[-irem]
      data_all[[dd]]$gev_year$lsl_max <- data_all[[dd]]$gev_year$lsl_max[-irem]
      data_all[[dd]]$gpd$year <- data_all[[dd]]$gpd$year[-irem]
      data_all[[dd]]$gpd$counts <- data_all[[dd]]$gpd$counts[-irem]
      data_all[[dd]]$gpd$excesses <- data_all[[dd]]$gpd$excesses[-irem]
      data_all[[dd]]$gpd$time_length <- data_all[[dd]]$gpd$time_length[-irem]
      data_all[[dd]]$gpd$time_length_all <- sum(data_all[[dd]]$gpd$time_length)
      data_all[[dd]]$gpd$counts_all <- sum(unlist(data_all[[dd]]$gpd$counts), na.rm=TRUE)
      data_all[[dd]]$gpd$excesses_all <- unlist(data_all[[dd]]$gpd$excesses)[!is.na(unlist(data_all[[dd]]$gpd$excesses))]
    } else if(data_all[[dd]]$gpd$year[1] > covariates[1,'year']) {
    # if begins after the forcing, clip the forcing
      irem_aux <- c(irem_aux, which(covariates[,'year'] < data_all[[dd]]$gpd$year[1]))
    }

    # if tide gauge record ends after auxiliary forcing, clip it
    if(max(data_all[[dd]]$gpd$year) > max(covariates[,'year'])) {
      irem <- which(data_all[[dd]]$gpd$year > max(covariates[,'year']))
      data_all[[dd]]$gev_year$year <- data_all[[dd]]$gev_year$year[-irem]
      data_all[[dd]]$gev_year$lsl_max <- data_all[[dd]]$gev_year$lsl_max[-irem]
      data_all[[dd]]$gpd$year <- data_all[[dd]]$gpd$year[-irem]
      data_all[[dd]]$gpd$counts <- data_all[[dd]]$gpd$counts[-irem]
      data_all[[dd]]$gpd$excesses <- data_all[[dd]]$gpd$excesses[-irem]
      data_all[[dd]]$gpd$time_length <- data_all[[dd]]$gpd$time_length[-irem]
      data_all[[dd]]$gpd$time_length_all <- sum(data_all[[dd]]$gpd$time_length)
      data_all[[dd]]$gpd$counts_all <- sum(unlist(data_all[[dd]]$gpd$counts), na.rm=TRUE)
      data_all[[dd]]$gpd$excesses_all <- unlist(data_all[[dd]]$gpd$excesses)[!is.na(unlist(data_all[[dd]]$gpd$excesses))]
    } else if(max(data_all[[dd]]$gpd$year) < max(covariates[,'year'])) {
    # if ends before forcing, clip the forcing
      irem_aux <- c(irem_aux, which(covariates[,'year'] > max(data_all[[dd]]$gpd$year)))
    }

    covariates_trimmed <- covariates
    if(length(irem_aux) > 0) {covariates_trimmed <- covariates_trimmed[-irem_aux,]}

    for (cc in names_covariates) {

      auxiliary_in <- covariates_trimmed[,cc]
      forc_max <- max(auxiliary_in)

      if(gpd.type==1) {
        auxiliary <- NULL
      } else {
        auxiliary <- auxiliary_in
      }

      out.deoptim <- DEoptim(neg_log_like_ppgpd, lower=gpd_models[[gpd.type]]$bound_lower, upper=gpd_models[[gpd.type]]$bound_upper,
                           DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                           parnames=gpd_models[[gpd.type]]$parnames, data_calib=data_all[[dd]]$gpd, auxiliary=auxiliary)
      deoptim.all[[cc]][[gpd.type]][dd,] <- out.deoptim$optim$bestmem
      colnames(deoptim.all[[cc]][[gpd.type]]) <- gpd_models[[gpd.type]]$parnames

    }
    tend <- proc.time()
    print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
  }
  tend0 <- proc.time()
  print(paste('... done. Took ',round(as.numeric(tend0-tbeg0)[3]/60,2),' minutes', sep=''))
}

print('...done.')

#
#===============================================================================
# check distributions, fit priors
#===============================================================================
#

print('fitting prior distributions to the MLE parameters...')

# fit gamma and normal priors
# -> centered at the medians
# -> with standard deviation equal to half the max-min range
#    (or do empirical sd? might underestimate though - take wider)

# assign which parameters have which priors
if(exists('gamma.priors')) {rm(list=c('gamma.priors','normal.priors','uniform.priors'))}
gamma.priors <- c('lambda','lambda0','sigma','sigma0')
normal.priors <- c('lambda1','sigma1','xi','xi0','xi1')
uniform.priors <- NULL

priors_normalgamma <- vector('list', length(names_covariates))
names(priors_normalgamma) <- names_covariates
for (cc in names_covariates) {
  priors_normalgamma[[cc]] <- vector('list', nmodel)
  for (model in 1:nmodel) {
    priors_normalgamma[[cc]][[model]] <- vector('list', length(gpd_models[[model]]$parnames))
    names(priors_normalgamma[[cc]][[model]]) <- gpd_models[[model]]$parnames
    for (par in gpd_models[[model]]$parnames) {
      priors_normalgamma[[cc]][[model]][[par]] <- vector('list', 3) # type, and 2 distribution parameters
      if(!is.na(match(par, uniform.priors))) {
         names(priors_normalgamma[[cc]][[model]][[par]]) <- c('type','shape','rate')
         priors_normalgamma[[cc]][[model]][[par]]$type <- 'uniform'
         priors_normalgamma[[cc]][[model]][[par]]$lower <- gpd_models[[model]]$bound_lower[match(par,gpd_models[[model]]$parnames)]
         priors_normalgamma[[cc]][[model]][[par]]$upper <- gpd_models[[model]]$bound_upper[match(par,gpd_models[[model]]$parnames)]
      } else if(!is.na(match(par, gamma.priors))) { # shape=alpha, rate=beta, mean=shape/rate, var=shape/rate^2
        names(priors_normalgamma[[cc]][[model]][[par]]) <- c('type','shape','rate')
        priors_normalgamma[[cc]][[model]][[par]]$type <- 'gamma'
        priors_normalgamma[[cc]][[model]][[par]]$rate <- mean(deoptim.all[[cc]][[model]][,par]) / var(deoptim.all[[cc]][[model]][,par])
        priors_normalgamma[[cc]][[model]][[par]]$shape <- mean(deoptim.all[[cc]][[model]][,par]) * priors_normalgamma[[cc]][[model]][[par]]$rate
      } else if(!is.na(match(par, normal.priors))) {
        names(priors_normalgamma[[cc]][[model]][[par]]) <- c('type','mean','sd')
        priors_normalgamma[[cc]][[model]][[par]]$type <- 'normal'
        priors_normalgamma[[cc]][[model]][[par]]$mean <- mean(deoptim.all[[cc]][[model]][,par])
        priors_normalgamma[[cc]][[model]][[par]]$sd   <- sd(deoptim.all[[cc]][[model]][,par])
      }
    }
  }
}

print('...done.')

#
#===============================================================================
# "fit" wide uniform priors (just using the bounds for the DE optim search)
#===============================================================================
#

print('fitting prior distributions to the uniform bounds for MLE search...')

# all parameters have uniform bounds, given by bound_lower_set and bound_upper_set
priors_uniform <- vector('list', length(names_covariates)); names(priors_uniform) <- names_covariates
for (cc in names_covariates) {
  priors_uniform[[cc]] <- vector('list', nmodel)
  for (model in 1:nmodel) {
    priors_uniform[[cc]][[model]] <- vector('list', length(gpd_models[[model]]$parnames)); names(priors_uniform[[cc]][[model]]) <- gpd_models[[model]]$parnames
    for (par in gpd_models[[model]]$parnames) {
      priors_uniform[[cc]][[model]][[par]] <- vector('list', 3) # type, and 2 distribution parameters
      names(priors_uniform[[cc]][[model]][[par]]) <- c('type','lower','upper'); priors_uniform[[cc]][[model]][[par]]$type <- 'uniform'
      priors_uniform[[cc]][[model]][[par]]$lower <- gpd_models[[model]]$bound_lower[match(par,gpd_models[[model]]$parnames)]
      priors_uniform[[cc]][[model]][[par]]$upper <- gpd_models[[model]]$bound_upper[match(par,gpd_models[[model]]$parnames)]
    }
  }
}

print('...done.')

#
#===============================================================================
# save priors
#===============================================================================
#

# rds -> save single object; the only one we need is 'priors'
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.priors.normalgamma <- paste(output.dir,'surge_priors_normalgamma_gpd_',today,'.rds', sep='')
filename.priors.uniform <- paste(output.dir,'surge_priors_uniform_gpd_',today,'.rds', sep='')
filename.mles <- paste(output.dir,'surge_MLEs_gpd_',today,'.rds', sep='')

print(paste('saving priors and DE optim output as .rds files to read and use later...',sep=''))

saveRDS(priors_normalgamma, file=filename.priors.normalgamma)
saveRDS(priors_uniform, file=filename.priors.uniform)
saveRDS(deoptim.all, file=filename.mles)
priors$gpd <- priors_normalgamma

print('...done.')



#===============================================================================
#===============================================================================
# GEV
#===============================================================================
#===============================================================================


#
#===============================================================================
# fit MLE GEV model parameters for each candidate model at each tide gauge
#===============================================================================
#

print('starting DE optimization for MLE GEV parameters for all stations in set...')

# need the likelihood function
source('likelihood_gev.R')

deoptim.all <- vector('list', length(names_covariates))
names(deoptim.all) <- names_covariates
for (cc in names_covariates) {
  deoptim.all[[cc]] <- vector('list', length(gev_models))
  for (i in 1:length(gev_models)) {
    deoptim.all[[cc]][[i]] <- mat.or.vec(length(data_all), length(gev_models[[i]]$parnames))
    rownames(deoptim.all[[cc]][[i]]) <- names(data_all)
    colnames(deoptim.all[[cc]][[i]]) <- gev_models[[i]]$parnames
  }
}

for (dd in 1:length(data_all)) {

  print(paste('starting to calculate MLE GEV parameters for tide gauge data set ',dd,' / ',length(data_all),sep=''))
  tbeg0 <- proc.time()
  data_all[[dd]]$deoptim <- vector('list', nmodel)

  for (gev.type in 1:length(gev_models)) {

    print(paste('  - starting DE optimization for model',gev.type,'...'))
    tbeg <- proc.time()

    irem_aux <- NULL

    # if tide gauge record starts before auxiliary forcing, clip it
    if(data_all[[dd]]$gev_year$year[1] < covariates[1,'year']) {
      irem <- which(data_all[[dd]]$gev_year$year < covariates[,'year'][1])
      data_all[[dd]]$gev_year$year <- data_all[[dd]]$gev_year$year[-irem]
      data_all[[dd]]$gev_year$lsl_max <- data_all[[dd]]$gev_year$lsl_max[-irem]
    } else if(data_all[[dd]]$gev_year$year[1] > covariates[1,'year']) {
      # if begins after the forcing, clip the forcing
      irem_aux <- c(irem_aux, which(covariates[,'year'] < data_all[[dd]]$gev_year$year[1]))
    }

    # if tide gauge record ends after auxiliary forcing, clip it
    if(max(data_all[[dd]]$gev_year$year) > max(covariates[,'year'])) {
      irem <- which(data_all[[dd]]$gev_year$year > max(covariates[,'year']))
      data_all[[dd]]$gev_year$year <- data_all[[dd]]$gev_year$year[-irem]
      data_all[[dd]]$gev_year$lsl_max <- data_all[[dd]]$gev_year$lsl_max[-irem]
    } else if(max(data_all[[dd]]$gev_year$year) < max(covariates[,'year'])) {
      # if ends before forcing, clip the forcing
      irem_aux <- c(irem_aux, which(covariates[,'year'] > max(data_all[[dd]]$gev_year$year)))
    }

    covariates_trimmed <- covariates
    if(length(irem_aux) > 0) {covariates_trimmed <- covariates_trimmed[-irem_aux,]}

    # now, match the years from the LSL_max record within the covariates' years
    idx_keep_aux <- match(data_all[[dd]]$gev_year$year, covariates_trimmed[,'year'])
    covariates_trimmed <- covariates_trimmed[idx_keep_aux,]

    for (cc in names_covariates) {

      auxiliary_in <- covariates_trimmed[,cc]
      forc_max <- max(auxiliary_in)

      if(gev.type==1) {
        auxiliary <- NULL
      } else {
        auxiliary <- auxiliary_in
      }

      out.deoptim <- DEoptim(neg_log_like_gev, lower=gev_models[[gev.type]]$bound_lower, upper=gev_models[[gev.type]]$bound_upper,
                           DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
                           parnames=gev_models[[gev.type]]$parnames, data_calib=data_all[[dd]]$gev_year$lsl_max, auxiliary=auxiliary)
      deoptim.all[[cc]][[gev.type]][dd,] <- out.deoptim$optim$bestmem
      colnames(deoptim.all[[cc]][[gev.type]]) <- gev_models[[gev.type]]$parnames

    }
    tend <- proc.time()
    print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
  }
  tend0 <- proc.time()
  print(paste('... done. Took ',round(as.numeric(tend0-tbeg0)[3]/60,2),' minutes', sep=''))
}

print('...done.')

#
#===============================================================================
# check distributions, fit priors
#===============================================================================
#

print('fitting prior distributions to the MLE parameters...')

# fit gamma and normal priors
# -> centered at the medians
# -> with standard deviation equal to half the max-min range
#    (or do empirical sd? might underestimate though - take wider)

# assign which parameters have which priors
if(exists('gamma.priors')) {rm(list=c('gamma.priors','normal.priors','uniform.priors'))}
gamma.priors <- c('mu','mu0','sigma','sigma0')
normal.priors <- c('mu1','sigma1','xi','xi0','xi1')
uniform.priors <- NULL

priors_normalgamma <- vector('list', length(names_covariates))
names(priors_normalgamma) <- names_covariates
for (cc in names_covariates) {
  priors_normalgamma[[cc]] <- vector('list', nmodel)
  for (model in 1:nmodel) {
    priors_normalgamma[[cc]][[model]] <- vector('list', length(gev_models[[model]]$parnames))
    names(priors_normalgamma[[cc]][[model]]) <- gev_models[[model]]$parnames
    for (par in gev_models[[model]]$parnames) {
      priors_normalgamma[[cc]][[model]][[par]] <- vector('list', 3) # type, and 2 distribution parameters
      if(!is.na(match(par, uniform.priors))) {
         names(priors_normalgamma[[cc]][[model]][[par]]) <- c('type','shape','rate')
         priors_normalgamma[[cc]][[model]][[par]]$type <- 'uniform'
         priors_normalgamma[[cc]][[model]][[par]]$lower <- gev_models[[model]]$bound_lower[match(par,gev_models[[model]]$parnames)]
         priors_normalgamma[[cc]][[model]][[par]]$upper <- gev_models[[model]]$bound_upper[match(par,gev_models[[model]]$parnames)]
      } else if(!is.na(match(par, gamma.priors))) { # shape=alpha, rate=beta, mean=shape/rate, var=shape/rate^2
        names(priors_normalgamma[[cc]][[model]][[par]]) <- c('type','shape','rate')
        priors_normalgamma[[cc]][[model]][[par]]$type <- 'gamma'
        priors_normalgamma[[cc]][[model]][[par]]$rate <- mean(deoptim.all[[cc]][[model]][,par]) / var(deoptim.all[[cc]][[model]][,par])
        priors_normalgamma[[cc]][[model]][[par]]$shape <- mean(deoptim.all[[cc]][[model]][,par]) * priors_normalgamma[[cc]][[model]][[par]]$rate
      } else if(!is.na(match(par, normal.priors))) {
        names(priors_normalgamma[[cc]][[model]][[par]]) <- c('type','mean','sd')
        priors_normalgamma[[cc]][[model]][[par]]$type <- 'normal'
        priors_normalgamma[[cc]][[model]][[par]]$mean <- mean(deoptim.all[[cc]][[model]][,par])
        priors_normalgamma[[cc]][[model]][[par]]$sd   <- sd(deoptim.all[[cc]][[model]][,par])
      }
    }
  }
}

print('...done.')

#
#===============================================================================
# "fit" wide uniform priors (just using the bounds for the DE optim search)
#===============================================================================
#

print('fitting prior distributions to the uniform bounds for MLE search...')

# all parameters have uniform bounds, given by bound_lower_set and bound_upper_set
priors_uniform <- vector('list', length(names_covariates)); names(priors_uniform) <- names_covariates
for (cc in names_covariates) {
  priors_uniform[[cc]] <- vector('list', nmodel)
  for (model in 1:nmodel) {
    priors_uniform[[cc]][[model]] <- vector('list', length(gev_models[[model]]$parnames)); names(priors_uniform[[cc]][[model]]) <- gev_models[[model]]$parnames
    for (par in gev_models[[model]]$parnames) {
      priors_uniform[[cc]][[model]][[par]] <- vector('list', 3) # type, and 2 distribution parameters
      names(priors_uniform[[cc]][[model]][[par]]) <- c('type','lower','upper'); priors_uniform[[cc]][[model]][[par]]$type <- 'uniform'
      priors_uniform[[cc]][[model]][[par]]$lower <- gev_models[[model]]$bound_lower[match(par,gev_models[[model]]$parnames)]
      priors_uniform[[cc]][[model]][[par]]$upper <- gev_models[[model]]$bound_upper[match(par,gev_models[[model]]$parnames)]
    }
  }
}

print('...done.')

#
#===============================================================================
# save priors
#===============================================================================
#

# rds -> save single object; the only one we need is 'priors'
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.priors.normalgamma <- paste(output.dir,'surge_priors_normalgamma_gev_',today,'.rds', sep='')
filename.priors.uniform <- paste(output.dir,'surge_priors_uniform_gev_',today,'.rds', sep='')
filename.mles <- paste(output.dir,'surge_MLEs_gev_',today,'.rds', sep='')

print(paste('saving priors and DE optim output as .rds files to read and use later...',sep=''))

saveRDS(priors_normalgamma, file=filename.priors.normalgamma)
saveRDS(priors_uniform, file=filename.priors.uniform)
saveRDS(deoptim.all, file=filename.mles)
priors$gev <- priors_normalgamma

print('...done.')

##==============================================================================
## End
##==============================================================================
