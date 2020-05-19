##==============================================================================
## write_excel_data_and_priors.R
##
## The main calibration and priors-fitting will create processed data and prior
## distribution RDS files. This routine will take those files and write Excel
## tables for them, so they are more easily read and shared. So, it is assumed
## that `process_data.R` and `fit_priors.R` will have been run before this.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

filename.processeddata.gev <- "../input_data/processeddata_gev_22Mar2020.rds"
filename.processeddata.gpd <- "../input_data/processeddata_gpd_23Mar2020.rds"
filename.covariates <- "../input_data/covariates_22Mar2020.rds"
filename.priors.normalgamma.gev <- "../surge_priors_normalgamma_gev_15May2020.rds"
filename.priors.normalgamma.gpd <- "../surge_priors_normalgamma_gpd_15May2020.rds"
filename.priors.uniform.gev <- "../surge_priors_uniform_gev_15May2020.rds"
filename.priors.uniform.gpd <- "../surge_priors_uniform_gpd_15May2020.rds"

todo

##==============================================================================
## End
##==============================================================================
