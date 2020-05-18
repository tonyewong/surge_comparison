##==============================================================================
## best_models.R
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================



##==============================================================================
## Maximum likelihood
##===================

## Read parameter results and covariate forcing
results <- vector('list', n_evm); names(results) <- names_evm
results$gev <- readRDS(filename_maxlike_gev)
results$gpd <- readRDS(filename_maxlike_gpd)
site_names <- names(results[[1]])
nsite <- length(site_names)

## Read covariate information and set up parameters for each model
data_calib <- vector('list', n_evm); names(data_calib) <- names_evm
data_calib$gev <- readRDS(filename_calibdata_gev)
data_calib$gpd <- readRDS(filename_calibdata_gpd)
covariates <- readRDS(filename_covariates)
names_covariates <- colnames(covariates)[2:5]
ncovar <- length(names_covariates)
source("trimmed_forcing.R")
source("parameter_setup_gev.R")
source("parameter_setup_gpd.R")
source("likelihood_gev.R")
source("likelihood_gpd.R")

##
## Calculate AIC, BIC and negative log-likelihood for each
## station (dim 1), for each covariate (dim 2), for each model (dim 3)
##

parameters <- nll <- aic <- bic <- vector('list', n_evm)
names(parameters) <- names(nll) <- names(aic) <- names(bic) <- names_evm

for (gg in names_evm) {
  if (gg=="gev")        {models <- gev_models
  } else if (gg=="gpd") {models <- gpd_models
  } else                {print("ERROR: unrecognized model type")
  }
  parameters[[gg]] <- vector('list', nmodel)
  nll[[gg]] <- aic[[gg]] <- bic[[gg]] <- array(dim=c(nsite,ncovar,nmodel), dimnames=list(site_names, names_covariates, 1:nmodel))
  for (mm in 1:nmodel) {
    parameters[[gg]][[mm]] <- vector("list", nsite)
    names(parameters[[gg]][[mm]]) <- site_names
    for (dd in site_names) {
      if (gg=="gev")        {ndata <- nrow(data_calib[[gg]][[dd]]); years <- data_calib[[gg]][[dd]][,"year"]
      } else if (gg=="gpd") {ndata <- data_calib[[gg]][[dd]]$counts_all; years <- data_calib[[gg]][[dd]]$year}
      parameters[[gg]][[mm]][[dd]] <- mat.or.vec(nr=ncovar, nc=length(models[[mm]]$parnames))
      rownames(parameters[[gg]][[mm]][[dd]]) <- names_covariates
      colnames(parameters[[gg]][[mm]][[dd]]) <- models[[mm]]$parnames
      for (cc in names_covariates) {
        covar_forc <- covariates[,cc]
        time_forc <- covariates[,"year"]
        parameters[[gg]][[mm]][[dd]][cc,] <- results[[gg]][[dd]][[cc]][[mm]]$optim$bestmem
        if (mm > 1) {auxiliary <- trimmed_forcing(years, time_forc, covar_forc)$forcing
        } else {auxiliary <- NULL}
        nll[[gg]][dd,cc,mm] <- results[[gg]][[dd]][[cc]][[mm]]$optim$bestval
        aic[[gg]][dd,cc,mm] <- 2*length(models[[mm]]$parnames) + 2*nll[[gg]][dd,cc,mm]
        bic[[gg]][dd,cc,mm] <- log(ndata)*length(models[[mm]]$parnames) + 2*nll[[gg]][dd,cc,mm]
      }
    }
  }
}

save(list=c("nll","aic","bic"), file="../output/comparison_metrics_maxlike.RData")
##==============================================================================



##==============================================================================
## Maximum a posteriori
##=====================

## Read parameter results and covariate forcing
results <- vector('list', n_evm); names(results) <- names_evm
results$gev <- readRDS(filename_maxpost_gev)
results$gpd <- readRDS(filename_maxpost_gpd)
site_names <- names(results[[1]])
nsite <- length(site_names)

##
## Calculate negative posterior score for each station (dim 1), for each
## covariate (dim 2), for each model (dim 3)
## Note: not computing AIC or BIC because the posterior score wraps up the
## probability of the parameters as a part of this (prior)
##

parameters <- nps <- vector('list', n_evm)
names(parameters) <- names(nps) <- names_evm

for (gg in names_evm) {
  if (gg=="gev")        {models <- gev_models
  } else if (gg=="gpd") {models <- gpd_models
  } else                {print("ERROR: unrecognized model type")
  }
  parameters[[gg]] <- vector('list', nmodel)
  nps[[gg]] <- array(dim=c(nsite,ncovar,nmodel), dimnames=list(site_names, names_covariates, 1:nmodel))
  for (mm in 1:nmodel) {
    parameters[[gg]][[mm]] <- vector("list", nsite)
    names(parameters[[gg]][[mm]]) <- site_names
    for (dd in site_names) {
      if (gg=="gev")        {ndata <- nrow(data_calib[[gg]][[dd]]); years <- data_calib[[gg]][[dd]][,"year"]
      } else if (gg=="gpd") {ndata <- data_calib[[gg]][[dd]]$counts_all; years <- data_calib[[gg]][[dd]]$year}
      parameters[[gg]][[mm]][[dd]] <- mat.or.vec(nr=ncovar, nc=length(models[[mm]]$parnames))
      rownames(parameters[[gg]][[mm]][[dd]]) <- names_covariates
      colnames(parameters[[gg]][[mm]][[dd]]) <- models[[mm]]$parnames
      for (cc in names_covariates) {
        covar_forc <- covariates[,cc]
        time_forc <- covariates[,"year"]
        parameters[[gg]][[mm]][[dd]][cc,] <- results[[gg]][[dd]][[cc]][[mm]]$optim$bestmem
        if (mm > 1) {auxiliary <- trimmed_forcing(years, time_forc, covar_forc)$forcing
        } else {auxiliary <- NULL}
        nps[[gg]][dd,cc,mm] <- results[[gg]][[dd]][[cc]][[mm]]$optim$bestval
      }
    }
  }
}

save(list=c("nps"), file="../output/comparison_metrics_maxpost.RData")
##==============================================================================



##==============================================================================
## End
##==============================================================================
