##==============================================================================
## likelihood_gev.R
##
## Likelihood function(s), priors, posterior for generalized extreme value
## distribution. Parameters may covary with some `auxiliary` time series.
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


##==============================================================================
## log(prior) for gev model
##==============================================================================

log_prior_gev <- function(parameters,
                          parnames,
                          priors
){
  lpri <- 0
  for (par in parnames) {
    parameter.value <- as.numeric(parameters[match(par,parnames)])
    if(priors[[par]]$type=='normal') {
      lpri <- lpri + dnorm(x=parameter.value, mean=priors[[par]]$mean, sd=priors[[par]]$sd, log=TRUE)
    } else if(priors[[par]]$type=='gamma') {
      lpri <- lpri + dgamma(x=parameter.value, shape=priors[[par]]$shape, rate=priors[[par]]$rate, log=TRUE)
    } else if(priors[[par]]$type=='uniform') {
      lpri <- lpri + dunif(x=parameter.value, min=priors[[par]]$lower, max=priors[[par]]$upper, log=TRUE)
    } else {print("ERROR: unknown priors type")}
  }
  return(lpri)
}
##==============================================================================


##==============================================================================
## log(likelihood) for gev model
##==============================================================================

log_like_gev <- function(parameters,     # set of GEV parameters (possibly nonstationary)
                         parnames,       # gev_models[[m]]$parnames
                         data_calib,     # processing_output[,"lsl_max"]
                         auxiliary=NULL  # time series of same length as data_calib
){
  if ("mu0" %in% parnames) {
    # location parameter nonstationary
    mu0 <- parameters[match('mu0',parnames)]
    mu1 <- parameters[match('mu1',parnames)]
    mu <- mu0 + mu1*auxiliary
  } else {
    mu <- parameters[match('mu',parnames)]
  }
  if ("sigma0" %in% parnames) {
    # scale parameter nonstationary
    sigma0 <- parameters[match('sigma0',parnames)]
    sigma1 <- parameters[match('sigma1',parnames)]
    sigma <- exp(sigma0 + sigma1*auxiliary)
  } else {
    sigma <- parameters[match('sigma',parnames)]
  }
  if ("xi0" %in% parnames) {
    # shape parameter nonstationary
    xi0 <- parameters[match('xi0',parnames)]
    xi1 <- parameters[match('xi1',parnames)]
    xi <- xi0 + xi1*auxiliary
  } else {
    xi <- parameters[match('xi',parnames)]
  }
  llik <- sum(devd(data_calib, loc=mu, scale=sigma, shape=xi, log=TRUE, type='GEV'))
  return(llik)
}
##==============================================================================


##==============================================================================
## negative log(likelihood) for gev model
##==============================================================================

neg_log_like_gev <- function(parameters,     # set of GEV parameters (possibly nonstationary)
                             parnames,       # gev_models[[m]]$parnames
                             data_calib,     # processing_output[,"lsl_max"]
                             auxiliary=NULL  # time series of same length as data_calib
){
  return(-log_like_gev(parameters,parnames,data_calib,auxiliary))
}
##==============================================================================


##==============================================================================
## log(posterior) for gev model
##==============================================================================

log_post_gev <- function(parameters, parnames, data_calib, priors, auxiliary){

  lpost <- 0
  llik <- 0
  lpri <- 0

  # calculate prior
  lpri <- log_prior_gev(parameters=parameters,
                        parnames=parnames,
                        priors=priors)

  if(is.finite(lpri)){
    # calculate likelihood (only if parameters pass the prior test)
    llik <- log_like_gev(parameters=parameters,
                         parnames=parnames,
                         data_calib=data_calib,
                         auxiliary=auxiliary)
  }

  lpost <- lpri + llik
  return(lpost)
}
##==============================================================================


##==============================================================================
## negative log(posterior) for gev model
##==============================================================================

neg_log_post_gev <- function(parameters, parnames, data_calib, priors, auxiliary) {
  return(-log_post_gev(parameters,parnames,data_calib,priors,auxiliary))
}
##==============================================================================


##==============================================================================
## End
##==============================================================================
