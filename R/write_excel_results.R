##==============================================================================
## write_excel_results.R
##
## The main calibration will create output results files from the optimization,
## including parameters for GEV and GPD that maximize the likelihood (or
## equivalently minimize the negative log-likelihood) and the posterior score.
## This routine will take those files and write Excel tables for them, so they
## are more easily read and shared. It is assumed that `calibration_driver.R`
## will have been run before this.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

filename.optimlike.gev <- "../output/optim_gev_like_22Mar2020.rds"
filename.optimpost.gev <- "../output/optim_gev_post_15May2020.rds"
filename.optimlike.gpd <- "../output/optim_gpd_like_23Mar2020.rds"
filename.optimpost.gpd <- "../output/optim_gpd_post_15May2020.rds"

todo

##==============================================================================
## End
##==============================================================================
