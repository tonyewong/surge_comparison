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

# write priors to CSV
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

# write processed calibration data to CSV

for (gg in evm_names) {
  for (dd in names(data_calib[[gg]])) {
    if (gg=="gev") {
      # simpler case
      filename.towrite <- paste("../csv/data/", gg, "/", dd, "_annmax.csv", sep="")
      write.csv(x=data_calib[[gg]][[dd]], file=filename.towrite, row.names=FALSE)
    } else if (gg=="gpd") {
      # not as simple case. need a bunch of info about the GPD
      filename.towrite <- paste("../csv/data/", gg, "/", dd, "_pot.csv", sep="")
      for (ii in 1:length(data_calib[[gg]][[dd]]$year)) {
        if (ii==1) {
          dat_to_write <- c(data_calib[[gg]][[dd]]$year[ii], data_calib[[gg]][[dd]]$counts[ii], data_calib[[gg]][[dd]]$time_length[ii], data_calib[[gg]][[dd]]$excesses[[ii]])
          dat_to_write <- t(data.frame(dat_to_write))
          col_names <- c("year","counts","time_length",rep("excesses", length(dat_to_write)-3))
          colnames(dat_to_write) <- col_names
          write.table(dat_to_write, file=filename.towrite, append=F, sep=',', row.names=F, col.names=T)
        } else {
          dat_to_write <- c(data_calib[[gg]][[dd]]$year[ii], data_calib[[gg]][[dd]]$counts[ii], data_calib[[gg]][[dd]]$time_length[ii], data_calib[[gg]][[dd]]$excesses[[ii]])
          dat_to_write <- t(data.frame(dat_to_write))
          write.table(dat_to_write, file=filename.towrite, append=T, sep=',', row.names=F, col.names=F)
        }
      }
    }
  }
}

# write the thresholds (in percent and in length units) for the POT exceedances, the declustering timescales

filename.towrite <- "../csv/data/gpd/thresholds_declustering.csv"
dat_to_write <- mat.or.vec(nr=length(data_calib$gpd), nc=3)
rownames(dat_to_write) <- names(data_calib$gpd)
colnames(dat_to_write) <- c("p.threshold", "threshold", "dt.decluster")
for (dd in names(data_calib$gpd)) {
  dat_to_write[dd,"p.threshold"] <- data_calib$gpd[[dd]]$p.threshold
  dat_to_write[dd,"threshold"] <- data_calib$gpd[[dd]]$threshold
  dat_to_write[dd,"dt.decluster"] <- data_calib$gpd[[dd]]$dt.decluster
}
write.csv(x=dat_to_write, file=filename.towrite)

##==============================================================================

if(FALSE){
# fixing p.threhsold / p.threshold typo
gg <- "gpd"
for (dd in names(data_calib[[gg]])) {
  data_calib[[gg]][[dd]]$p.threshold <- data_calib[[gg]][[dd]]$p.threhsold
  data_calib[[gg]][[dd]]$p.threhsold <- NULL
}
saveRDS(data_calib[[gg]], file=filename.processeddata[[gg]])
}

##==============================================================================
## End
##==============================================================================
