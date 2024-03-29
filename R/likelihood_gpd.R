##==============================================================================
## likelihood_gpd.R
##
## likelihood function(s), priors, posterior
## for generalized pareto/poisson process distribution
##
## The idea is that a GPD governs how excesses above a given threshold are
## distributed, but this is conditioned on the probability that we observe an
## excess. These exceedance observations are governed by a Poisson process, where
## the rate parameter gives (1/) the expected number of exceedances per year (or
## other time unit desired; here I use a year).
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
## pdf and cdf for pp-gpd (that can handle arrays)
## note that this assumes you aren't a jerk and send in different length
## arrays. so... don't do that.
##==============================================================================

# use eq. 11 of Martins and Stedinger (2001), for each year independently:
# llike <- log([likelihood of seeing exactly m exceedances in time interval dt with rate lambda[t]])
#           + log([joint GPD density for the m exceedances in time interval dt])
# ^^^ do this for all of the intervals

# need to send in for each interval:
# 1. length of intervals
# 2. number of exceedances each interval (assumed to be [level-threshold])
# 3. all parameters (rate, scale, shape) as time series (annual values/each block)
# 4. values of the exceedance levels (level-threshold)


# assumes that x comes in as x-threshold
# also doesn't address the case where some of the shape=0 and others do not
gpd_pdf <- function(x, scale, shape, log=FALSE){
  p <- rep(NA,length(x))
  if(all(shape==0)) {
    if(log) {
      p <- -x/scale - log(scale)
    } else {
      p <- (1/scale) * exp(-x/scale)
    }
  } else {
    if(log) {
      p <- -(1+1/shape)*log(1+shape*x/scale) - log(scale)
      if(any(shape < 0)) {
        i1 <- which(shape < 0); i2 <- which(x > (-scale/shape)); i3 <- intersect(i1,i2)
        p[i3] <- -Inf
      } else if(any(shape > 0)) {
        i1 <- which(shape > 0); i2 <- which(x < (-scale/shape)); i3 <- intersect(i1,i2)
        p[i3] <- -Inf
      }
    } else {
      p <- (1/scale) * (1+shape*x/scale)^(-(1+1/shape))
      if(any(shape < 0)) {
        i1 <- which(shape < 0); i2 <- which(x > (loc-scale/shape)); i3 <- intersect(i1,i2)
        p[i3] <- 0
      } else if(any(shape > 0)) {
        i1 <- which(shape > 0); i2 <- which(x < (loc-scale/shape)); i3 <- intersect(i1,i2)
        p[i3] <- 0
      }
    }
  }
  return(p)
}

gpd_cdf <- function(q, scale, shape){
  p <- rep(NA,length(q))
  if(all(shape==0)) {
    p <- 1-exp(-q/scale)
  } else {
    p <- 1-((1+shape*q/scale)^(-1/shape))
    if(any(shape < 0)) {
      i1 <- which(shape < 0); i2 <- which(q > (-scale/shape)); i3 <- intersect(i1,i2)
      # if shape < 0, i3 is all the places where q > theoretical upper bound
      # -> there ought to be 100% probability mass below here
      p[i3] <- 1
    } else if(any(shape > 0)) {
      i1 <- which(shape > 0); i2 <- which(q < (-scale/shape)); i3 <- intersect(i1,i2)
      # if shape > 0, i3 is all the places where q < theoretical upper bound
      # -> there ought to be 100% probability mass above here
      p[i3] <- 0
    }
  }
  return(p)
}

# h is effective height (mm), lambda is poisson process rate (lambda*365.25 is
# expected number of exceedances/year), scale and shape are as in GPD, threshold
# is the GPD threshold, nmax is the maximum number of events/year we consider,
# time.length is the number of days in the time interval considered (365.25 for
# a year).

ppgpd_overtop <- function(h, lambda, sigma, xi, threshold, nmax, time.length){
  res <- 0
  #cdf_gpd <- pevd(q=h-threshold, threshold=0, scale=sigma, shape=xi, type='GP')
  cdf_gpd <- gpd_cdf(q=(h-threshold), scale=sigma, shape=xi)
  pdf_pp <- dpois(x=1:nmax, lambda=(lambda*time.length))
  for (i in 1:nmax) {
    res <- res + pdf_pp[i] * (1-(cdf_gpd^i))
  }
  return(res)
}
##==============================================================================



##==============================================================================
## project PP-GPD parameters
##==============================================================================

project_ppgpd <- function(parameters,
                          parnames,
                          auxiliary
){
  parameters_project <- mat.or.vec(length(auxiliary), 3)
  colnames(parameters_project) <- c('lambda','sigma','xi')
  n.param <- length(parnames)
  if ("lambda0" %in% parnames) {
    # Poisson process rate parameter nonstationary
    lambda0 <- parameters[match('lambda0',parnames)]
    lambda1 <- parameters[match('lambda1',parnames)]
    lambda <- lambda0 + lambda1*auxiliary
  } else {
    lambda <- rep(parameters[match('lambda',parnames)], nbins)
  }
  if ("sigma0" %in% parnames) {
    # scale parameter nonstationary
    sigma0 <- parameters[match('sigma0',parnames)]
    sigma1 <- parameters[match('sigma1',parnames)]
    sigma <- exp(sigma0 + sigma1*auxiliary)
  } else {
    sigma <- rep(parameters[match('sigma',parnames)], nbins)
  }
  if ("xi0" %in% parnames) {
    # shape parameter nonstationary
    xi0 <- parameters[match('xi0',parnames)]
    xi1 <- parameters[match('xi1',parnames)]
    xi <- xi0 + xi1*auxiliary
  } else {
    xi <- rep(parameters[match('xi',parnames)], nbins)
  }

  parameters_project[,'lambda'] <- lambda
  parameters_project[,'sigma'] <- sigma
  parameters_project[,'xi'] <- xi

  return(parameters_project)
}
##==============================================================================



##==============================================================================
## log(prior) for PP-GPD model
##==============================================================================

log_prior_ppgpd <- function(parameters,
                            parnames,
                            priors,
                            model,
                            auxiliary=NULL
){
  lpri <- 0

  for (par in parnames) {
    parameter.value <- as.numeric(parameters[match(par,parnames)])
    if(priors[[model]][[par]]$type=='normal') {
      lpri <- lpri + dnorm(x=parameter.value, mean=priors[[model]][[par]]$mean, sd=priors[[model]][[par]]$sd, log=TRUE)
    } else if(priors[[model]][[par]]$type=='gamma') {
      lpri <- lpri + dgamma(x=parameter.value, shape=priors[[model]][[par]]$shape, rate=priors[[model]][[par]]$rate, log=TRUE)
    } else if(priors[[model]][[par]]$type=='uniform') {
      lpri <- lpri + dunif(x=parameter.value, min=priors[[model]][[par]]$lower, max=priors[[model]][[par]]$upper, log=TRUE)
    }
  }

  return(lpri)
}
##==============================================================================



##==============================================================================
## -log(likelihood) for PP-GPD model
##==============================================================================

neg_log_like_ppgpd <- function(parameters,
                               parnames,
                               data_calib,
                               auxiliary=NULL
){
  nll <- -1 * log_like_ppgpd(parameters, parnames, data_calib, auxiliary)
  return(nll)
}
##==============================================================================



##==============================================================================
## log(likelihood) for PP-GPD model
##==============================================================================

log_like_ppgpd <- function(parameters,
                           parnames,
                           data_calib,
                           auxiliary=NULL
){
  llik <- 0
  nbins <- length(data_calib$counts)
  #print(parameters)
  n.param <- length(parnames)
  if ("lambda0" %in% parnames) {
    # Poisson process rate parameter nonstationary
    lambda0 <- parameters[match('lambda0',parnames)]
    lambda1 <- parameters[match('lambda1',parnames)]
    lambda <- lambda0 + lambda1*auxiliary
  } else {
    lambda <- rep(parameters[match('lambda',parnames)], nbins)
  }
  if ("sigma0" %in% parnames) {
    # scale parameter nonstationary
    sigma0 <- parameters[match('sigma0',parnames)]
    sigma1 <- parameters[match('sigma1',parnames)]
    sigma <- exp(sigma0 + sigma1*auxiliary)
  } else {
    sigma <- rep(parameters[match('sigma',parnames)], nbins)
  }
  if ("xi0" %in% parnames) {
    # shape parameter nonstationary
    xi0 <- parameters[match('xi0',parnames)]
    xi1 <- parameters[match('xi1',parnames)]
    xi <- xi0 + xi1*auxiliary
  } else {
    xi <- rep(parameters[match('xi',parnames)], nbins)
  }

  # check extra conditions that would otherwise 'break' the likelihood function
  if( any(lambda < 0) | any(sigma < 0) | any(is.na(lambda)) | any(is.na(sigma)) ) {llik <- -Inf}
  else {
    llik.bin <- dpois(x=data_calib$counts, lambda=(lambda*data_calib$time_length), log=TRUE)
    hits <- which(data_calib$counts!=0)   # gets indices of bins (years) with exceedances
    for (b in hits) {
      llik.bin[b] <- llik.bin[b] +
                     sum(devd(data_calib$excesses[[b]]-data_calib$threshold, threshold=0, scale=sigma[b], shape=xi[b], log=TRUE, type='GP'))
      #if( (xi[b] < 0) & any( (data_calib$excesses[[b]]-data_calib$threshold) < (-sigma[b]/xi[b]) )) {llik.bin[b] <- -Inf; break;}
    }
    llik <- sum(llik.bin)
  }

  # constraint on lambda0, lambda1 and forc_max
  if( exists('lambda0') & exists('lambda1') ) {
    forc_max <- max(auxiliary)
    if( lambda1[1] < (-lambda0[1]/forc_max) ) {llik <- -Inf}
  }

  return(llik)
}
##==============================================================================



##==============================================================================
## log(post) for PP-GPD model
##==============================================================================

log_post_ppgpd <- function(parameters,
                           parnames,
                           data_calib,
                           priors,
                           model,
                           auxiliary
){
  lpost <- 0
  llik <- 0
  lpri <- 0

  # calculate prior
  lpri <- log_prior_ppgpd(parameters=parameters,
                          parnames=parnames,
                          priors=priors,
                          model=model,
                          auxiliary=auxiliary)

  if(is.finite(lpri)){
    # calculate likelihood (only if parameters pass the prior test)
    llik <- log_like_ppgpd(parameters=parameters,
                           parnames=parnames,
                           data_calib=data_calib,
                           auxiliary=auxiliary)
  }

  lpost <- lpri + llik
  return(lpost)
}
##==============================================================================


##==============================================================================
## negative log(posterior) for gpd model
##==============================================================================

neg_log_post_ppgpd <- function(parameters, parnames, data_calib, priors, model, auxiliary) {
  return(-log_post_ppgpd(parameters,parnames,data_calib,priors,model,auxiliary))
}
##==============================================================================


##==============================================================================
## End
##==============================================================================
