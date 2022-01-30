##==============================================================================
## write_csv_parameterss.R
##
## The main calibration will create output results files from the optimization,
## including parameters for GEV and GPD that maximize the likelihood (or
## equivalently minimize the negative log-likelihood) and the posterior score.
## This routine will take those files and write CSV tables for them, so they
## are more easily read and shared. It is assumed that `calibration_driver.R`
## will have been run before this.
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
parameters <- filename.optim <- vector('list', 2)
names(parameters) <- names(filename.optim) <- c("nll","nps")
for (metric in names(parameters)) {
  parameters[[metric]] <- filename.optim[[metric]] <- vector('list', 2)
  names(parameters[[metric]]) <- names(filename.optim[[metric]]) <- c("gev","gpd")
}
source("parameter_setup_gev.R")
source("parameter_setup_gpd.R")

# set the file names that you used
filename.optim$nll$gev <- "../output/optim_gev_like_28Jan2022.rds"
filename.optim$nll$gpd <- "../output/optim_gpd_like_28Jan2022.rds"
filename.optim$nps$gev <- "../output/optim_gev_post_28Jan2022.rds"
filename.optim$nps$gpd <- "../output/optim_gpd_post_28Jan2022.rds"

# read the calibration output
for (metric in names(parameters)) {
  for (gg in evm_names) {
    parameters[[metric]][[gg]] <- readRDS(filename.optim[[metric]][[gg]])
  }
}

# write parameters to files
for (metric in names(parameters)) {
  for (gg in evm_names) {
    if (gg=="gev") models <- gev_models
    if (gg=="gpd") models <- gpd_models
    for (cc in names(parameters[[metric]][[gg]][[1]])) {
      for (mm in 1:length(models)) {
        parameters_towrite <- mat.or.vec(nr=length(parameters[[metric]][[gg]]), nc=length(parameters[[metric]][[gg]][[1]][[1]][[mm]]$optim$bestmem))
        rownames(parameters_towrite) <- names(parameters[[metric]][[gg]])
        colnames(parameters_towrite) <- models[[mm]]$parnames
        for (dd in names(parameters[[metric]][[gg]])) {
          parameters_towrite[dd,] <- parameters[[metric]][[gg]][[dd]][[cc]][[mm]]$optim$bestmem
        }
        filename.towrite <- paste("../csv/parameters/", metric, "/", gg, "/", cc, "/model", mm, ".csv", sep="")
        write.csv(x=parameters_towrite, file=filename.towrite)
      }
    }
  }
}

##==============================================================================
## End
##==============================================================================
