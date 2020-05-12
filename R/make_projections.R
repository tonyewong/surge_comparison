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
rl <- rp <- vector('list', n_evm); names(rl) <- names(rp) <- names_evm
for (gg in names_evm) {
  rl[[gg]] <- rp[[gg]] <- vector('list', nsite); names(rl[[gg]]) <- names(rp[[gg]]) <- site_names
  for (dd in 1:nsite) {
    rl[[gg]][[dd]] <- rp[[gg]][[dd]] <- vector("list", ncovar)
    names(rl[[gg]][[dd]]) <- names(rp[[gg]][[dd]]) <- names_covariates
    for (cc in names_covariates) {
      rl[[gg]][[dd]][[cc]] <- rp[[gg]][[dd]][[cc]] <- array(dim=c(nrp, nmodel, nyear), dimnames=list(return_periods, 1:nmodel, years_proj))
    }
  }
}

for (gg in names_evm) {
  for (dd in 1:nsite) {
    print(paste("Making ",gg," projections for site ",site_names[dd]," (",dd,"/",length(data_calib[[gg]]),")...", sep=""))
    for (cc in names_covariates) {
      covariate_proj <- covariates[match(years_proj, covariates[,'year']), cc]
      # this is pretty clunky. can replace later but let's be honest, you probably won't and anyone who's reading this probably also is like "yeah, me neither" -TW
      if (gg=="gev") {
        for (mm in 1:nmodel) {
          params <- parameters[[gg]][[mm]][[dd]][cc,]
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
              rl[[gg]][[dd]][[cc]][rr,mm,yy] <- qevd(1-1/return_periods[rr], loc=mu, scale=sigma, shape=xi)
              if (yy==1) {
                rp[[gg]][[dd]][[cc]][rr,mm,yy] <- 1/(1-pevd(rl[[gg]][[dd]][[cc]][rr,mm,yy], loc=mu, scale=sigma, shape=xi))
              } else {
                rp[[gg]][[dd]][[cc]][rr,mm,yy] <- 1/(1-pevd(rl[[gg]][[dd]][[cc]][rr,mm,1], loc=mu, scale=sigma, shape=xi))
              }
            }
          }
        }
      } else if (gg=="gpd") {
        for (mm in 1:nmodel) {
          params <- parameters[[gg]][[mm]][[dd]][cc,]
          parnames <- gpd_models[[mm]]$parnames
          for (yy in 1:nyear) {
            if ("lambda0" %in% parnames) {
              # rate parameter nonstationary
              lambda0 <- params[match('lambda0',parnames)]
              lambda1 <- params[match('lambda1',parnames)]
              lambda <- lambda0 + lambda1*covariate_proj[yy]
            } else {
              lambda <- params[match('lambda',parnames)]
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
              rl[[gg]][[dd]][[cc]][rr,mm,yy] <- rlevd(return_periods[rr], scale=sigma,
                                                      shape=xi, threshold=data_calib[[gg]][[dd]]$threshold,
                                                      type='GP', npy=365.25, rate=lambda)
              if (yy==1) {
                rp[[gg]][[dd]][[cc]][rr,mm,yy] <- 1/ppgpd_overtop(rl[[gg]][[dd]][[cc]][rr,mm,yy], lambda, sigma, xi,
                                                                  data_calib[[gg]][[dd]]$threshold, nmax=182.6, time.length=365.25)
              } else {
                rp[[gg]][[dd]][[cc]][rr,mm,yy] <- 1/ppgpd_overtop(rl[[gg]][[dd]][[cc]][rr,mm,1], lambda, sigma, xi,
                                                                  data_calib[[gg]][[dd]]$threshold, nmax=182.6, time.length=365.25)
              }
            }
          }
        }
      }
    }
  }
}

## Save return levels object
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.returnlevels <- paste('../output/returnlevels_',today,'.RData', sep='')
print(paste("Saving return levels data (`rl`) to file",filename.returnlevels))
save(list=c('rl', 'rp', 'site_names', 'covariates', 'return_periods', 'years_proj',
            'nsite', 'nmodel', 'nrp', 'nyear', 'gev_models', 'gpd_models', 'data_calib'), file=filename.returnlevels)

##==============================================================================
## End
##==============================================================================
