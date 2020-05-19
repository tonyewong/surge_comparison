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

parameters_maxlike <- nll <- aic <- bic <- vector('list', n_evm)
names(parameters_maxlike) <- names(nll) <- names(aic) <- names(bic) <- names_evm

for (gg in names_evm) {
  if (gg=="gev")        {models <- gev_models
  } else if (gg=="gpd") {models <- gpd_models
  } else                {print("ERROR: unrecognized model type")
  }
  parameters_maxlike[[gg]] <- vector('list', nmodel)
  nll[[gg]] <- aic[[gg]] <- bic[[gg]] <- array(dim=c(nsite,ncovar,nmodel), dimnames=list(site_names, names_covariates, 1:nmodel))
  for (mm in 1:nmodel) {
    parameters_maxlike[[gg]][[mm]] <- vector("list", nsite)
    names(parameters_maxlike[[gg]][[mm]]) <- site_names
    for (dd in site_names) {
      if (gg=="gev")        {ndata <- nrow(data_calib[[gg]][[dd]]); years <- data_calib[[gg]][[dd]][,"year"]
      } else if (gg=="gpd") {ndata <- data_calib[[gg]][[dd]]$counts_all; years <- data_calib[[gg]][[dd]]$year}
      parameters_maxlike[[gg]][[mm]][[dd]] <- mat.or.vec(nr=ncovar, nc=length(models[[mm]]$parnames))
      rownames(parameters_maxlike[[gg]][[mm]][[dd]]) <- names_covariates
      colnames(parameters_maxlike[[gg]][[mm]][[dd]]) <- models[[mm]]$parnames
      for (cc in names_covariates) {
        covar_forc <- covariates[,cc]
        time_forc <- covariates[,"year"]
        parameters_maxlike[[gg]][[mm]][[dd]][cc,] <- results[[gg]][[dd]][[cc]][[mm]]$optim$bestmem
        if (mm > 1) {auxiliary <- trimmed_forcing(years, time_forc, covar_forc)$forcing
        } else {auxiliary <- NULL}
        nll[[gg]][dd,cc,mm] <- results[[gg]][[dd]][[cc]][[mm]]$optim$bestval
        aic[[gg]][dd,cc,mm] <- 2*length(models[[mm]]$parnames) + 2*nll[[gg]][dd,cc,mm]
        bic[[gg]][dd,cc,mm] <- log(ndata)*length(models[[mm]]$parnames) + 2*nll[[gg]][dd,cc,mm]
      }
    }
  }
}
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

parameters_maxpost <- nps <- vector('list', n_evm)
names(parameters_maxpost) <- names(nps) <- names_evm

for (gg in names_evm) {
  if (gg=="gev")        {models <- gev_models
  } else if (gg=="gpd") {models <- gpd_models
  } else                {print("ERROR: unrecognized model type")
  }
  parameters_maxpost[[gg]] <- vector('list', nmodel)
  nps[[gg]] <- array(dim=c(nsite,ncovar,nmodel), dimnames=list(site_names, names_covariates, 1:nmodel))
  for (mm in 1:nmodel) {
    parameters_maxpost[[gg]][[mm]] <- vector("list", nsite)
    names(parameters_maxpost[[gg]][[mm]]) <- site_names
    for (dd in site_names) {
      if (gg=="gev")        {ndata <- nrow(data_calib[[gg]][[dd]]); years <- data_calib[[gg]][[dd]][,"year"]
      } else if (gg=="gpd") {ndata <- data_calib[[gg]][[dd]]$counts_all; years <- data_calib[[gg]][[dd]]$year}
      parameters_maxpost[[gg]][[mm]][[dd]] <- mat.or.vec(nr=ncovar, nc=length(models[[mm]]$parnames))
      rownames(parameters_maxpost[[gg]][[mm]][[dd]]) <- names_covariates
      colnames(parameters_maxpost[[gg]][[mm]][[dd]]) <- models[[mm]]$parnames
      for (cc in names_covariates) {
        covar_forc <- covariates[,cc]
        time_forc <- covariates[,"year"]
        parameters_maxpost[[gg]][[mm]][[dd]][cc,] <- results[[gg]][[dd]][[cc]][[mm]]$optim$bestmem
        if (mm > 1) {auxiliary <- trimmed_forcing(years, time_forc, covar_forc)$forcing
        } else {auxiliary <- NULL}
        nps[[gg]][dd,cc,mm] <- results[[gg]][[dd]][[cc]][[mm]]$optim$bestval
      }
    }
  }
}
##==============================================================================



##==============================================================================
## Save results to an RData file
##==============================

save(list=c("nll","aic","bic","nps"), file="../output/comparison_metrics.RData")
##==============================================================================



##==============================================================================
## End
##==============================================================================
