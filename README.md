# surge_comparison

1. process_data.R
  1. process_gev.R - yields `processeddata_gev_[date].rds`
  1. get_timeseries_covariates.R - yields `covariates_[date].rds`
1. calibration_driver.R
  1. trimmed_forcing.R
  1. likelihood_gev.R
  1. parameter_setup_gev.R
  1. yields `optim_[covariate name]-gev-[date].rds`, output from maximum likelihood optimization for GEV parameters in each of the 8 potentially nonstationary models.
