##==============================================================================
## make_projections.R
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

## What return periods and years do you want?
return_periods <- c(2,5,10,20,50,100,200,500)
years_proj <- seq(from=2020, to=2099)
nyear <- length(years_proj)
nrp <- length(return_periods)


## Make projections using each of the models

# want projections for each site, for each model, for a variety of return periods, for a variety of years, for each simulation in the ensemble
# Note: set up like this, can reference the year we want as string in 3rd element: rl[[dd]][,,"2050"], e.g.
rl <- vector('list', n_evm); names(rl) <- names_evm
for (gg in names_evm) {
  rl[[gg]] <- vector('list', nsite); names(rl[[gg]]) <- site_names
  for (dd in 1:nsite) {
    rl[[gg]][[dd]] <- vector("list", ncovar)
    names(rl[[gg]][[dd]]) <- names_covariates
    for (cc in names_covariates) {
      rl[[gg]][[dd]][[cc]] <- array(dim=c(nrp, nmodel, nyear), dimnames=list(return_periods, 1:nmodel, years_proj))
    }
  }
}

TODO - modify to include gpd models

for (gg in names_evm) {
  for (dd in 1:nsite) {
    print(paste("Making projections for site ",site_names[dd]," (",dd,"/",length(data_calib),")...", sep=""))
    for (cc in names_covariates) {
      covariate_proj <- covariates[match(years_proj, covariates[,'year']), cc]
      # this is pretty clunky. can replace later but let's be honest, you probably won't and anyone who's reading this probably also is like "yeah, me neither" -TW
      if (gg=="gev") {
        for (mm in 1:nmodel) {
          params <- parameters[[mm]][[dd]][cc,]
          parnames <- gev_models[[mm]]$parnames
          for (yy in 1:nyear) {
            if ("mu0" %in% parnames) {
              # location parameter nonstationary
              mu0 <- params[match('mu0',parnames)]
              mu1 <- params[match('mu1',parnames)]
              mu <- mu0 + mu1*covariate_proj[yy]
            } else {
              mu <- params[match('mu',parnames)]
            }
            if ("sigma0" %in% parnames) {
              # scale parameter nonstationary
              sigma0 <- params[match('sigma0',parnames)]
              sigma1 <- params[match('sigma1',parnames)]
              sigma <- exp(sigma0 + sigma1*covariate_proj[yy])
            } else {
              sigma <- params[match('sigma',parnames)]
            }
            if ("xi0" %in% parnames) {
              # shape parameter nonstationary
              xi0 <- params[match('xi0',parnames)]
              xi1 <- params[match('xi1',parnames)]
              xi <- xi0 + xi1*covariate_proj[yy]
            } else {
              xi <- params[match('xi',parnames)]
            }
            for (rr in 1:nrp) {
              rl[[dd]][[cc]][rr,mm,yy] <- qevd(1-1/return_periods[rr], loc=mu, scale=sigma, shape=xi)
            }
          }
        }
      } else if (gg=="gpd") {

      }
    }
  }
}

## Save return levels object
filename.returnlevels <- paste('../output/returnlevels_',sim_id,'.RData', sep='')
save(list=c('rl', 'site_names', 'covariates', 'return_periods', 'years_proj',
            'nsite', 'nmodel', 'nrp', 'nyear', 'sim_id'), file=filename.returnlevels)

##==============================================================================
## End
##==============================================================================
