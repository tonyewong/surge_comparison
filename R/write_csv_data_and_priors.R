##==============================================================================
## write_csv_data_and_priors.R
##
## The main calibration and priors-fitting will create processed data and prior
## distribution RDS files. This routine will take those files and write CSV
## tables for them, so they are more easily read and shared. So, it is assumed
## that `process_data.R` and `fit_priors.R` will have been run before this.
##
## Questions? Tony Wong (aewsma@rit.edu)
##==============================================================================

# initialize
filename.priors <- vector('list', 2); names(filename.priors) <- c("normalgamma","uniform")
filename.processeddata <- filename.priors$normalgamma <- filename.priors$uniform <- vector('list', 2)
names(filename.processeddata) <- names(filename.priors.normalgamma) <- names(filename.priors.uniform) <- evm_names <- c("gev", "gpd")

# set the file names that you used
filename.processeddata$gev <- "../input_data/processeddata_gev_22Mar2020.rds"
filename.processeddata$gpd <- "../input_data/processeddata_gpd_23Mar2020.rds"
filename.priors$normalgamma$gev <- "../input_data/surge_priors_normalgamma_gev_15May2020.rds"
filename.priors$normalgamma$gpd <- "../input_data/surge_priors_normalgamma_gpd_15May2020.rds"
filename.priors$uniform$gev <- "../input_data/surge_priors_uniform_gev_15May2020.rds"
filename.priors$uniform$gpd <- "../input_data/surge_priors_uniform_gpd_15May2020.rds"

# read the covariates
filename.covariates <- "../input_data/covariates_22Mar2020.rds"
covariates <- readRDS(filename.covariates)

# read the processed calibration data
data_calib <- vector('list', 2); names(data_calib) <- evm_names
for (gg in evm_names) {
  data_calib[[gg]] <- readRDS(filename.processeddata[[gg]])
}

# read the prior distributions
priors <- vector('list', length(filename.priors)); names(priors) <- names(filename.priors)
priors$normalgamma <- priors$uniform <- vector('list', 2)
names(priors$normalgamma) <- names(priors$uniform) <- evm_names
for (distr in names(filename.priors)) {
  for (gg in evm_names) {
    priors[[distr]][[gg]] <- readRDS(filename.priors[[distr]][[gg]])
  }
}

# write covariates to CSV
write.csv(x=covariates, file="../csv/covariates.csv", row.names=FALSE)

# write priors to files
for (distr in names(filename.priors)) {
  for (gg in evm_names) {
    for (cc in names(priors[[distr]][[gg]])) {
      for (mm in 1:length(priors[[distr]][[gg]][[cc]])) {
        filename.towrite <- paste("../csv/priors/", distr, "/", gg, "/", cc, "/model", mm, ".csv", sep="")
        write.csv(x=unlist(priors[[distr]][[gg]][[cc]][[mm]]), file=filename.towrite)
      }
    }
  }
}

##==============================================================================
## End
##==============================================================================
