# surge_comparison

## Workflow

1. process_data.R
  1. process_gev.R - yields `processeddata_gev_[date].rds`
  1. process_gpd.R - yields `processeddata_gpd_[date].rds`; note that this processing takes a while (hours) because of the moving average to account for mean sea level rise
  1. get_timeseries_covariates.R - yields `covariates_[date].rds`
1. calibration_driver.R
  1. trimmed_forcing.R
  1. likelihood_gev.R
  1. likelihood_gpd.R
  1. parameter_setup_gev.R
  1. parameter_setup_gpd.R
  1. yields `optim_[covariate name]-gev-[date].rds`, output from maximum likelihood optimization for GEV parameters in each of the 8 potentially nonstationary models.


## Input data

### Tide gauge stations

* Only take stations with at least 15 years of available data


