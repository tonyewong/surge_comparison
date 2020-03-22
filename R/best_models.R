##==============================================================================
## best_models.R
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

## Read parameter results and covariate forcing
results <- readRDS(filename_optim)
site_names <- names(results)
nsite <- length(site_names)

## Read covariate information and set up parameters for each model
data_calib <- readRDS(filename_calibdata)
covariates <- readRDS(filename_covariates)
names_covariates <- colnames(covariates)[2:5]
ncovar <- length(names_covariates)
source("trimmed_forcing.R")
source("parameter_setup_gev.R")
source("likelihood_gev.R")

##
## Calculate AIC, BIC and negative log-likelihood for each
## station (dim 1), for each covariate (dim 2), for each model (dim 3)
##

nll <- aic <- bic <- array(dim=c(nsite,ncovar,nmodel), dimnames=list(site_names, names_covariates, 1:nmodel))
parameters <- vector('list', nmodel)

for (mm in 1:nmodel) {
  parameters[[mm]] <- vector("list", nsite)
  names(parameters[[mm]]) <- site_names
  for (dd in site_names) {
    parameters[[mm]][[dd]] <- mat.or.vec(nr=ncovar, nc=length(gev_models[[mm]]$parnames))
    rownames(parameters[[mm]][[dd]]) <- names_covariates
    colnames(parameters[[mm]][[dd]]) <- gev_models[[mm]]$parnames
    for (cc in names_covariates) {
      covar_forc <- covariates[,cc]
      time_forc <- covariates[,"year"]
      parameters[[mm]][[dd]][cc,] <- results[[dd]][[cc]][[mm]]$optim$bestmem
      if (mm > 1) {auxiliary <- trimmed_forcing(data_calib[[dd]][,"year"], time_forc, covar_forc)$forcing
      } else {auxiliary <- NULL}
      nll[dd,cc,mm] <- neg_log_like_gev(parameters[[mm]][[dd]][cc,], parnames=gev_models[[mm]]$parnames, data_calib=data_calib[[dd]][,"lsl_max"], auxiliary=auxiliary)
      aic[dd,cc,mm] <- 2*length(gev_models[[mm]]$parnames) + 2*nll[dd,cc,mm]
      bic[dd,cc,mm] <- log(nrow(data_calib[[dd]]))*length(gev_models[[mm]]$parnames) + 2*nll[dd,cc,mm]
    }
  }
}

write.csv(rbind(rep("NLL",nmodel), nll, rep("AIC",nmodel), aic, rep("BIC",nmodel), bic), file="../output/comparison_metrics.csv")



## Save return levels object
#filename.returnlevels <- paste('../output/returnlevels_',sim_id,'.RData', sep='')
#save(list=c('rl', 'bma_weights', 'site_names', 'covariates', 'return_periods',
#            'years_proj', 'nsite', 'nmodel', 'nrp', 'nyear', 'covariate_proj',
#            'sim_id'), file=filename.returnlevels)

##==============================================================================
## End
##==============================================================================
