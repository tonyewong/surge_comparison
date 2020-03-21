##==============================================================================
## info_content.R
##
## Tony Wong (aewsma@rit.edu)
##==============================================================================


se <- kl <- vector('list', nsite); names(se) <- names(kl) <- site_names
for (dd in 1:nsite) {
  se[[dd]] <- kl[[dd]] <- array(dim=c(nmodel, nyear), dimnames=list(1:nmodel, years_proj))
}

print("Making projections of return levels.")
x <- seq(from=-5000, to=25000, by=1)
for (dd in 1:nsite) {
  print(paste("Starting on date set", site_names[dd],"..."))
  # to compute KL divergence, need a baseline. use stationary distribution
  fx0 <- devd(x, loc=parameters[[1]][dd,"mu"], scale=parameters[[1]][dd,"sigma"], shape=parameters[[1]][dd,"xi"], type="GEV")
  for (mm in 1:nmodel) {
    params <- parameters[[mm]][dd,]
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
      fx <- devd(x, loc=mu, scale=sigma, shape=xi, type="GEV")
      # compute differential entropy, analogous to shannon/kl but continuous:
      se[[dd]][mm,yy] <- sintegral(x,-fx*log(fx))$value
      kl[[dd]][mm,yy] <- sintegral(x,fx0*log(fx0/fx))$value
    }
  }
}



##==============================================================================
## End
##==============================================================================
