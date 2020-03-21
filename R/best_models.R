##==============================================================================
## best_models.R
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

## Set up. Make sure consistent with the
distr <- "gev"
datestamp <- "21Mar2020"
sim_id <- paste(distr,datestamp, sep="-")
filename_optim <- paste("../output/optim_",sim_id,".rds", sep="")
filename_covariates <- "../input_data/covariates_21Mar2020.rds"
filename_calibdata  <- "../input_data/processeddata_gev_21Mar2020.rds"

## Read parameter results and covariate forcing
results <- readRDS(filename_optim)
site_names <- names(results)
nsite <- length(site_names)

## Read covariate information and set up parameters for each model
data_calib <- readRDS(filename_calibdata)
covariates <- readRDS(filename_covariates)
source("trimmed_forcing.R")
source("parameter_setup_gev.R")
source("likelihood_gev.R")

##
## Calculate AIC, BIC and negative log-likelihood for each model, for each station
##

nll <- aic <- bic <- mat.or.vec(nr=nsite, nc=nmodel)
rownames(nll) <- rownames(aic) <- rownames(bic) <- site_names
parameters <- vector('list', nmodel)

# TODO -- loop over covariates too - construct a 3D array
covar_forc <- covariates[,covar_name]
time_forc <- covariates[,"year"]

for (mm in 1:nmodel) {
  parameters[[mm]] <- mat.or.vec(nr=nsite, nc=length(gev_models[[mm]]$parnames))
  rownames(parameters[[mm]]) <- site_names
  colnames(parameters[[mm]]) <- gev_models[[mm]]$parnames
  for (dd in site_names) {
    parameters[[mm]][dd,] <- results[[dd]][[mm]]$optim$bestmem
    if (mm > 1) {auxiliary <- trimmed_forcing(data_calib[[dd]][,"year"], time_forc, covar_forc)$forcing
    } else {auxiliary <- NULL}
    nll[dd,mm] <- neg_log_like_gev(parameters[[mm]][dd,], parnames=gev_models[[mm]]$parnames, data_calib=data_calib[[dd]][,"lsl_max"], auxiliary=auxiliary)
    aic[dd,mm] <- 2*length(gev_models[[mm]]$parnames) + 2*nll[dd,mm]
    bic[dd,mm] <- log(nrow(data_calib[[dd]]))*length(gev_models[[mm]]$parnames) + 2*nll[dd,mm]
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
